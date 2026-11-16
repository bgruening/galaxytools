# pynxtools_dataconverter

Converts one raw file and its ELN to a NeXus file with
[pynxtools](https://github.com/FAIRmat-NFDI/pynxtools) and a reader plugin, chosen in the `reader`
select. Today: `spm` from [pynxtools-spm](https://github.com/FAIRmat-NFDI/pynxtools-spm), in the
public pynxtools-spm container. A new reader is one more case of that select, with its definitions.
A file that cannot be converted gives a status JSON with the reason; the job still succeeds.

## Raw file header

`raw_file_header` writes a raw file's header as text for the ELN tools. It reads the file with
[gwyddionpy](https://github.com/FAIRmat-NFDI/gwyddionPy) (many SPM formats, also binary ones) in the
same container; files gwyddionpy cannot read give their leading text.

## Tests

| Files | Checks |
|---|---|
| `const_dos_No14_003.sxm`, `eln_data.yaml` | A Nanonis STM image converts to NXstm |
| `broken.sxm`, `raw_without_extension` | Recorded failures; the job stays green |
| `eln_no_technique.yaml` | `auto` without a definition is a recorded failure |
| `const_dos_No14_003.sxm` (header tool) | gwyddionpy reads the Nanonis header |
| `sts_header_cut.dat` (header tool) | Nanonis spectrum: the text before `[DATA]`; cut from a pynxtools-spm test file (Apache-2.0) |
| `broken.sxm` (header tool) | Unknown format: its leading text |
