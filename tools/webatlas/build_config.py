#!/usr/bin/env python3
"""
build_config.py
====================================
Generates a Vitessce View config
"""

from __future__ import annotations
import typing as T
from collections import defaultdict
import os
import fire
import json
import logging
from itertools import chain, cycle
import regex
from vitessce import (
    VitessceConfig,
    DataType as dt,
    FileType as ft,
    CoordinationType as ct,
    Component as cm,
)
from constants.constants import (
    DATA_TYPES,
    DEFAULT_OPTIONS,
    DEFAULT_LAYOUTS,
    COMPONENTS_COORDINATION_TYPES,
    COMPONENTS_DATA_TYPES,
)


def build_options(
    file_type: str,
    file_path: str,
    file_options: dict[str, T.Any],
    check_exist: bool = False,
) -> T.Any:
    """Function that creates the View config's options for non-image files

    Args:
        file_type (str): Type of file supported by Vitessce.
        file_path (str): Path to file.
        file_options (dict[str, T.Any]): Dictionary defining the options.
        check_exist (bool, optional): Whether to check the given path to confirm the file exists.
            Defaults to False.

    Returns:
        T.Any: Options dictionary for View config file
    """
    options = None
    dts = set([])

    # @TODO: Change yaml input keys to 'locations', 'embeddings', 'labels'
    if file_type == ft.ANNDATA_ZARR:
        options = {}
        if "spatial" in file_options and "xy" in file_options["spatial"]:
            values_path = os.path.join(file_path, file_options["spatial"]["xy"])
            if not (check_exist and not os.path.exists(values_path)):
                options["obsLocations"] = {"path": file_options["spatial"]["xy"]}
                dts.add(dt.OBS_LOCATIONS)

        if "mappings" in file_options:
            for k, v in file_options["mappings"].items():
                values_path = os.path.join(file_path, k)
                if check_exist and not os.path.exists(values_path):
                    continue
                embedding_type = k.split("/")[-1].upper()
                v = v if v is not None and v != "null" else [0, 1]
                embedding = {"path": k, "embeddingType": embedding_type, "dims": v}
                options.setdefault("obsEmbedding", []).append(embedding)
                dts.add(dt.OBS_EMBEDDING)

        if "factors" in file_options:
            for factor in file_options["factors"]:
                values_path = os.path.join(file_path, factor)
                if check_exist and not os.path.exists(values_path):
                    continue
                label_type = factor.split("/")[-1].capitalize()
                label = {"path": factor, "obsLabelsType": label_type}
                options.setdefault("obsLabels", []).append(label)
                dts.add(dt.OBS_LABELS)

        if "sets" in file_options:
            for obs_set in file_options["sets"]:
                if type(obs_set) != dict:
                    obs_set = {"name": obs_set}
                if check_exist:
                    if not os.path.exists(os.path.join(file_path, obs_set["name"])) or (
                        "score" in obs_set
                        and not os.path.exists(
                            os.path.join(file_path, obs_set["score"])
                        )
                    ):
                        continue
                obs_set_options = {
                    "name": "".join(
                        w.capitalize()
                        for w in obs_set["name"].split("/")[-1].split("_")
                    ),
                    "path": obs_set["name"],
                }
                if "score" in obs_set:
                    obs_set_options["name"] += " with scores"
                    obs_set_options["scorePath"] = obs_set["score"]
                options.setdefault("obsSets", []).append(obs_set_options)
                dts.add(dt.OBS_SETS)

        if "matrix" in file_options:
            options["obsFeatureMatrix"] = {"path": file_options["matrix"]}
            dts.add(dt.OBS_FEATURE_MATRIX)

    return options, dts


def build_raster_options(
    images: dict[str, list[dict[str, T.Any]]], url: str
) -> dict[str, T.Any]:
    """Function that creates the View config's options for image files

    Args:
        images (dict[str, list[dict[str, T.Any]]], optional): Dictionary containing for each image type key (raw and label)
            a list of dictionaries (one per image of that type) with the corresponding path and metadata for that image.
            Defaults to {}.
        url (str): URL to prepend to each file in the config file.
            The URL to the local or remote server that will serve the files

    Returns:
        dict[str, T.Any]: Options dictionary for View config file
    """
    raster_options = {"renderLayers": [], "schemaVersion": "0.0.2", "images": []}
    for img_type in images.keys():  # raw, label
        for img in images[img_type]:
            image_name = img.get("name")
            channel_names = (
                img["md"]["channel_names"]
                if "channel_names" in img["md"] and len(img["md"]["channel_names"])
                else (
                    ["Labels"]
                    if img_type == "label"
                    else [f"Channel {x}" for x in range(int(img["md"]["C"]))]
                )
            )
            isBitmask = img_type == "label"
            raster_options["renderLayers"].append(image_name)
            raster_options["images"].append(
                {
                    "name": image_name,
                    "url": f"{url}{img['path']}",
                    "type": "zarr",
                    "metadata": {
                        "isBitmask": isBitmask,
                        "dimensions": [
                            {"field": "t", "type": "quantitative", "values": None},
                            {
                                "field": "channel",
                                "type": "nominal",
                                "values": channel_names,
                            },
                            {"field": "y", "type": "quantitative", "values": None},
                            {"field": "x", "type": "quantitative", "values": None},
                        ],
                        "isPyramid": True,
                        "transform": {"translate": {"y": 0, "x": 0}, "scale": 1},
                    },
                }
            )
    return raster_options


def write_json(
    project: str = "",
    dataset: str = "",
    file_paths: list[str] = [],
    images: dict[str, list[dict[str, T.Any]]] = {},
    url: str = "",
    options: dict[str, T.Any] = None,
    layout: str = "minimal",
    custom_layout: str = None,
    title: str = "",
    description: str = "",
    config_filename_suffix: str = "config.json",
    outdir: str = "./",
) -> None:
    """This function writes a Vitessce View config JSON file

    Args:
        project (str, optional): Project name. Defaults to "".
        dataset (str, optional): Dataset name. Defaults to "".
        file_paths (list[str], optional): Paths to files that will be included in the config file. Defaults to [].
        images (dict[str, list[dict[str, T.Any]]], optional): Dictionary containing for each image type key (raw and label)
            a list of dictionaries (one per image of that type) with the corresponding path and metadata for that image.
            Defaults to {}.
        url (str, optional): URL to prepend to each file in the config file.
            The URL to the local or remote server that will serve the files.
            Defaults to "".
        options (dict[str, T.Any], optional): Dictionary with Vitessce config file `options`. Defaults to None.
        layout (str, optional): Type of predefined layout to use. Defaults to "minimal".
        custom_layout (str, optional): String defining a Vitessce layout following its alternative syntax.
            https://vitessce.github.io/vitessce-python/api_config.html#vitessce.config.VitessceConfig.layout
            https://github.com/vitessce/vitessce-python/blob/1e100e4f3f6b2389a899552dffe90716ffafc6d5/vitessce/config.py#L855
            Defaults to None.
        title (str, optional): Data title to show in the visualization. Defaults to "".
        config_filename_suffix (str, optional): Config filename suffix. Defaults to "config.json".
        outdir (str, optional): Directory in which the config file will be written to. Defaults to "./".

    Raises:
        SystemExit: If no valid files have been input
        SystemExit: If the layout has an error that can not be fixed
    """
    if isinstance(file_paths, str):
        file_paths = json.loads(file_paths)

    if isinstance(images, str):
        images = json.loads(images)

    if isinstance(options, str):
        options = json.loads(options)

    has_files = False

    config = VitessceConfig(
        "1.0.15",
        name=str(title) if len(title) else str(project),
        description=description,
    )
    config_dataset = config.add_dataset(str(dataset), str(dataset))

    coordination_types = defaultdict(lambda: cycle(iter([])))

    dts = set([])

    if images.keys() and any([len(images[k]) for k in images.keys()]):
        has_files = True
        config_dataset.add_file(
            ft.RASTER_JSON, options=build_raster_options(images, url)
        )
        dts.add(dt.RASTER)

    for file_path in file_paths:
        file_type = ft.ANNDATA_ZARR

        has_files = True
        if file_type in DEFAULT_OPTIONS:
            file_options, file_dts = build_options(
                file_type, file_path, options or DEFAULT_OPTIONS[file_type]
            )
            # Set a coordination scope for any 'obsEmbedding'
            if "obsEmbedding" in file_options:
                coordination_types[ct.EMBEDDING_TYPE] = cycle(
                    chain(
                        coordination_types[ct.EMBEDDING_TYPE],
                        [
                            config.set_coordination_value(
                                ct.EMBEDDING_TYPE.value,
                                k["embeddingType"],
                                k["embeddingType"],
                            )
                            for k in file_options["obsEmbedding"]
                        ],
                    )
                )
        else:
            file_options, file_dts = None, []

        if file_options and len(file_options):
            config_dataset.add_file(
                file_type,
                url=f"{url}{file_path}",
                options=file_options,
            )
            dts.update(file_dts)
            break

    if not has_files:
        raise SystemExit("No files to add to config file")

    # Get layout components/views
    # Set layout with alternative syntax:
    # https://vitessce.github.io/vitessce-python/api_config.html#vitessce.config.VitessceConfig.layout
    # https://github.com/vitessce/vitessce-python/blob/1e100e4f3f6b2389a899552dffe90716ffafc6d5/vitessce/config.py#L855
    config_layout = (
        custom_layout
        if custom_layout and len(custom_layout)
        else DEFAULT_LAYOUTS[layout]
    )
    views, views_indx = [], []
    for m in regex.finditer("[a-zA-Z]+", config_layout):
        try:
            component = cm(m.group())
        except ValueError:
            logging.warning("{} is not a valid component".format(m.group()))
            views_indx.append((m.start(), m.end(), ""))
            continue

        # Remove component from layout if its required data type is not present
        if component in COMPONENTS_DATA_TYPES and COMPONENTS_DATA_TYPES[
            component
        ].isdisjoint(dts):
            logging.warning(
                "No file of data type required by component {}".format(component)
            )
            views_indx.append((m.start(), m.end(), ""))
            continue

        view = config.add_view(component, dataset=config_dataset)

        if component in COMPONENTS_COORDINATION_TYPES:
            has_ctype = True
            for coordination_type in COMPONENTS_COORDINATION_TYPES[component]:
                # Cycle through coordination scopes with the required coordination type
                try:
                    view.use_coordination(next(coordination_types[coordination_type]))
                except StopIteration:
                    logging.warning(
                        "No coordination scope with coordination type {} for component {}".format(
                            coordination_type, component
                        )
                    )
                    has_ctype = False
                    break
            if not has_ctype:
                views_indx.append((m.start(), m.end(), ""))
                continue

        views.append(view)
        views_indx.append((m.start(), m.end(), "views[{}]".format(len(views) - 1)))

    # Replace components names in layout string with the corresponding view object
    for start, end, view in sorted(views_indx, key=lambda x: x[0], reverse=True):
        config_layout = config_layout[:start] + view + config_layout[end:]

    # Remove remaining | and / operators and empty parenthesis
    empty_par = r"\((?>[^\(\)a-zA-Z]|(?R))*\)"
    unbalanced_op = r"(?<=\()([\/\|]+)|([\/\|]+)(?=\))|(?<=^)([\/\|]+)|([\/\|]+)(?=$)"
    repeated_op = r"(?<=[\/\|])([\/\|]+)"
    config_layout = regex.sub(empty_par, "", config_layout)
    config_layout = regex.sub(unbalanced_op, "", config_layout)
    config_layout = regex.sub(repeated_op, "", config_layout)

    # Concatenate views
    try:
        exec("config.layout({})".format(config_layout))
    except Exception as e:
        logging.error(e)
        raise SystemExit("Error building config layout")

    # Make sure views' grid coordinates are integers
    config_json = config.to_dict()
    for l in config_json["layout"]:
        for k in ("x", "y", "w", "h"):
            l[k] = round(l[k])

    if outdir and not os.path.isdir(outdir):
        os.mkdir(outdir)
    with open(
        os.path.join(outdir or "", f"{project}-{dataset}-{config_filename_suffix}"), "w"
    ) as out_file:
        json.dump(config_json, out_file, indent=2)


if __name__ == "__main__":
    fire.Fire(write_json)
