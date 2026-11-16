# pynxtools_dataconverter

Converts one raw file and its ELN to a NeXus file with
[pynxtools](https://github.com/FAIRmat-NFDI/pynxtools) and the `spm` reader of
[pynxtools-spm](https://github.com/FAIRmat-NFDI/pynxtools-spm), in the public pynxtools-spm container.
A file that cannot be converted gives a status JSON with the reason; the job still succeeds.

## Tests

| Files | Checks |
|---|---|
| `const_dos_No14_003.sxm`, `eln_data.yaml` | A Nanonis STM image converts to NXstm |
| `broken.sxm`, `raw_without_extension` | Recorded failures; the job stays green |
| `eln_no_technique.yaml` | `auto` without a definition is a recorded failure |
