# pynxtools_spm

The tools that need [pynxtools-spm](https://github.com/FAIRmat-NFDI/pynxtools-spm). They share one
environment and one version, the pynxtools-spm version. Another pynxtools plugin gets its own folder
like this one.

## Converter

`pynxtools_spm_converter` converts one raw file and its ELN to a NeXus file with
[pynxtools](https://github.com/FAIRmat-NFDI/pynxtools) and the `spm` reader.
A file that cannot be converted gives a status JSON with the reason; the job still succeeds.

## Raw file header

`raw_file_header` writes a raw file's header (metadata) as JSON for the ELN tools. It reads the file with
[gwyddionpy](https://github.com/FAIRmat-NFDI/gwyddionPy) (many SPM formats, also binary ones) in the
same environment; files gwyddionpy cannot read give their leading text.

## ELN files

`pynxtools_spm_eln_files` writes the plugin's ELN rules and instructions (`src/pynxtools_spm/galaxy/`,
from pynxtools-spm 0.5.1) for the ELN tools in `tools/nexus_eln`, once per workflow run. Optional
override files are used instead, e.g. fixes to test before they go into pynxtools-spm.

## Tests

| Files | Checks |
|---|---|
| `const_dos_No14_003.sxm`, `eln_data.yaml` | A Nanonis STM image converts to NXstm |
| `broken.sxm`, `raw_without_extension` | Recorded failures; the job stays green |
| `eln_no_technique.yaml` | `auto` without a definition is a recorded failure |
| `const_dos_No14_003.sxm` (header tool) | gwyddionpy reads the Nanonis header |
| `sts_header_cut.dat` (header tool) | Nanonis spectrum: the text before `[DATA]`; cut from a pynxtools-spm test file (Apache-2.0) |
| `broken.sxm` (header tool) | Unknown format: its leading text |
| none (ELN files tool) | The package's rules and instructions |
| `eln_rules_override.yaml`, `eln_instructions_override.md` (ELN files tool) | Overrides are passed on; one override alone |
| `eln_rules_no_definitions.yaml` (ELN files tool) | Rules without `definitions` fail |
