#!/usr/bin/env python
"""Peak finding for an already-converted XPS NeXus file.

Reads the default NXdata group's signal/axes (energy vs. intensity) out of an
existing NeXus/HDF5 file, runs ``scipy.signal.find_peaks`` on the intensity, and
appends one NXfit/NXpeak group per detected peak carrying the peak's energy
(``data/position``) and intensity (``data/intensity``).

This is a standalone port of ``pynxtools_xps.peak_finding`` that writes the same
NeXus structure using plain h5py, so the Galaxy tool depends only on scipy and
h5py (both on conda-forge). Upstream lives on the unreleased
``eosc-galaxy-peak-finding`` branch of pynxtools-xps and additionally requires
pynxtools, which is not packaged for conda. See:
https://github.com/FAIRmat-NFDI/nomad-eosc-galaxy-actions/tree/main/examples/galaxy_demo

The NXfit written here is deliberately partial: ``scipy.signal.find_peaks``
locates peaks but does not fit a background or a parametric peak shape, so
NXfit's own required ``data`` (fit input/output) and ``backgroundBACKGROUND``
groups are left out rather than filled with invented values. ``pynx validate``
will flag both as missing -- that reflects this step's actual scope (peak
finding, not peak fitting), not a bug.
"""

import argparse
import json
import logging
import numbers
import shutil
import sys
from pathlib import Path

import h5py
from scipy.signal import find_peaks

__version__ = "0.1.0"

FIT_LABEL = "Peaks found via scipy.signal.find_peaks"

logger = logging.getLogger("xps_peak_finding")


def decode_if_string(value):
    """Return *value* as ``str`` if it is a (possibly HDF5 byte) string."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, str):
        return value
    return value


def nx_class(node):
    """Return the ``NX_class`` attribute of *node* as ``str``, or None."""
    return decode_if_string(node.attrs.get("NX_class"))


def find_default_nxentry(root):
    """Return the default NXentry under *root*, else the first NXentry found."""
    default = decode_if_string(root.attrs.get("default"))
    if default and default in root:
        candidate = root[default]
        if isinstance(candidate, h5py.Group) and nx_class(candidate) == "NXentry":
            return candidate

    for key in root:
        child = root[key]
        if isinstance(child, h5py.Group) and nx_class(child) == "NXentry":
            return child
    return None


def _follow_default_chain(group):
    """Follow ``@default`` from *group* until it stops pointing at a group."""
    current = group
    while True:
        default = decode_if_string(current.attrs.get("default"))
        if not default or default not in current:
            return current if current is not group else None
        target = current[default]
        if not isinstance(target, h5py.Group):
            return current if current is not group else None
        current = target


def find_default_nxdata(nxentry):
    """Return the default plottable NXdata group under *nxentry*.

    Prefers the NIAC2014 (v3) ``@default`` chain, which must land on a group
    carrying ``@signal``; falls back to the first NXdata child (v2/v1 layouts).
    """
    nxdata = _follow_default_chain(nxentry)
    if nxdata is not None and decode_if_string(nxdata.attrs.get("signal")):
        return nxdata

    for key in nxentry:
        child = nxentry[key]
        if isinstance(child, h5py.Group) and nx_class(child) == "NXdata":
            return child
    return None


def inspect_nxdata(group):
    """Return ``(signal, first_axis)`` datasets for the NXdata *group*.

    Implements the v3 convention (group-level ``@signal``/``@axes``) first, then
    falls back to v2 (a field carrying its own ``signal="1"`` and a
    colon-separated ``@axes``). Returns ``(None, None)`` if neither applies.
    """
    # v3 (NIAC2014): the group names its signal and axes.
    signal_name = decode_if_string(group.attrs.get("signal"))
    if signal_name and signal_name in group:
        signal = group[signal_name]
        axes = group.attrs.get("axes")
        if isinstance(axes, (str, bytes)):
            axis_names = [decode_if_string(axes)]
        elif axes is not None:
            axis_names = [decode_if_string(axis) for axis in axes]
        else:
            axis_names = []
        for name in axis_names:
            if name in group:
                return signal, group[name]

    # v2: a field flags itself as the signal and lists its axes.
    for key in group:
        dataset = group[key]
        if not isinstance(dataset, h5py.Dataset):
            continue
        if decode_if_string(dataset.attrs.get("signal")) != "1":
            continue
        axes = decode_if_string(dataset.attrs.get("axes"))
        if isinstance(axes, str):
            for name in axes.replace(",", ":").split(":"):
                if name and name in group:
                    return dataset, group[name]

    # v1: axis fields carry an integer @axis attribute instead.
    for key in group:
        dataset = group[key]
        if not isinstance(dataset, h5py.Dataset):
            continue
        if decode_if_string(dataset.attrs.get("signal")) != "1":
            continue
        for axis_key in group:
            axis = group[axis_key]
            if isinstance(axis, h5py.Dataset) and isinstance(
                axis.attrs.get("axis"), numbers.Integral
            ):
                return dataset, axis

    return None, None


def _write_scalar(group, name, value, units):
    """Write a scalar dataset, attaching ``@units`` when the source had them."""
    dataset = group.create_dataset(name, data=value)
    if units is not None:
        dataset.attrs["units"] = decode_if_string(units)


def find_xps_peaks(input_path, output_path=None, **find_peaks_kwargs):
    """Find peaks in an XPS spectrum and append them to the NeXus file.

    Args:
        input_path: path to an existing NeXus file with a resolvable default
            NXdata group (energy vs. intensity).
        output_path: where to write the result. Defaults to ``input_path``
            (in place). If different, ``input_path`` is copied first and left
            untouched.
        **find_peaks_kwargs: forwarded to ``scipy.signal.find_peaks`` (e.g.
            ``height``, ``prominence``, ``distance``).

    Returns:
        One dict per detected peak, with ``position`` (energy) and ``intensity``.
    """
    input_path = Path(input_path)
    output_path = Path(output_path) if output_path is not None else input_path
    if output_path != input_path:
        shutil.copyfile(input_path, output_path)

    with h5py.File(output_path, "r+") as nexus_file:
        entry = find_default_nxentry(nexus_file)
        if entry is None:
            raise ValueError(f"No NXentry found in {output_path}")

        nxdata_group = find_default_nxdata(entry)
        if nxdata_group is None:
            raise ValueError(f"No default NXdata group found in {output_path}")

        signal, axis = inspect_nxdata(nxdata_group)
        if signal is None or axis is None:
            raise ValueError(
                f"NXdata group '{nxdata_group.name}' in {output_path} has no "
                "resolvable signal/axes."
            )

        intensity = signal[()]
        energy = axis[()]
        position_units = axis.attrs.get("units")
        intensity_units = signal.attrs.get("units")

        peak_indices, _properties = find_peaks(intensity, **find_peaks_kwargs)

        # NXfit is a recommended child group directly under ENTRY (a sibling of
        # the spectrum's own NXdata group), not nested under it -- see
        # NXxps.nxdl.xml. Drop any previous run so re-running is idempotent.
        if "fit" in entry:
            del entry["fit"]
        fit_group = entry.create_group("fit")
        fit_group.attrs["NX_class"] = "NXfit"
        fit_group.create_dataset("label", data=FIT_LABEL)

        peaks = []
        for i, index in enumerate(peak_indices):
            position = float(energy[index])
            peak_intensity = float(intensity[index])
            peaks.append({"position": position, "intensity": peak_intensity})

            peak_id = f"peak_{i}"
            peak_group = fit_group.create_group(peak_id)
            peak_group.attrs["NX_class"] = "NXpeak"
            peak_group.create_dataset("label", data=peak_id)
            data_group = peak_group.create_group("data")
            data_group.attrs["NX_class"] = "NXdata"
            _write_scalar(data_group, "position", position, position_units)
            _write_scalar(data_group, "intensity", peak_intensity, intensity_units)

    return peaks


def write_tsv(peaks, path):
    """Write *peaks* as a tab-separated table with a header row."""
    with open(path, "w") as handle:
        handle.write("peak_id\tposition\tintensity\n")
        for i, peak in enumerate(peaks):
            handle.write(f"peak_{i}\t{peak['position']}\t{peak['intensity']}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Find peaks in an XPS NeXus spectrum and append them as NXfit/NXpeak."
    )
    parser.add_argument(
        "input_path", help="Path to an already-converted XPS NeXus file."
    )
    parser.add_argument(
        "-o", "--output-path", default=None, help="Output path (defaults to in-place)."
    )
    parser.add_argument("--height", type=float, default=None)
    parser.add_argument("--prominence", type=float, default=None)
    parser.add_argument("--distance", type=float, default=None)
    parser.add_argument(
        "--peaks-tsv", default=None, help="Also write the peak table as TSV."
    )
    parser.add_argument(
        "--peaks-json", default=None, help="Also write the peak table as JSON."
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )

    kwargs = {
        key: value
        for key, value in {
            "height": args.height,
            "prominence": args.prominence,
            "distance": args.distance,
        }.items()
        if value is not None
    }

    peaks = find_xps_peaks(args.input_path, args.output_path, **kwargs)
    logger.info("Found %d peaks in %s", len(peaks), args.input_path)

    if args.peaks_tsv:
        write_tsv(peaks, args.peaks_tsv)
    if args.peaks_json:
        with open(args.peaks_json, "w") as handle:
            json.dump(peaks, handle, indent=2)


if __name__ == "__main__":
    main()
