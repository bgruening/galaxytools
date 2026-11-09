# Zenodo harvest for SPM data

Goal: turn openly licensed scanning probe microscopy (SPM) datasets on [Zenodo](https://zenodo.org)
into [NeXus](https://www.nexusformat.org) files and [NOMAD](https://nomad-lab.eu) projects that link
back to the Zenodo DOI.

The `zenodo_harvest` Action of
[nomad-eosc-galaxy-actions](https://github.com/FAIRmat-NFDI/nomad-eosc-galaxy-actions) decides what to
harvest, starts this workflow and creates the NOMAD projects from its outputs. The workflow only
computes; it writes nothing to NOMAD.

```
work list ─► download + unzip ─► ELN builder ─► NeXus converter ─► report builder
                                  (per raw file)  (per raw file)
```

The tools are in [`tools/zenodo_nomad_harvest`](../../tools/zenodo_nomad_harvest/).

## Inputs and outputs

| Name | Content |
|---|---|
| `work list` (input, step 0) | `work_list.json` from the Action |
| `prefetched files` (optional input) | Files already in Galaxy; element identifier = Zenodo file name |
| `manifest` | Every archive and file, with its role and what happened to it |
| `results` | Conversion status per raw file |
| `projects` | One NOMAD project zip per Zenodo record |
| `logs`, `summary` | A log per record, and an HTML summary of the run |
| `eln files`, `eln reports`, `nexus files`, `conversion statuses`, `converter logs` | Outputs of the middle steps |

The Action reads the outputs by these labels, so keep them when you edit the workflow.

## Files

| File | Use |
|---|---|
| `zenodo-harvest-spm.ga` | Import in Galaxy (*Workflows → Import*, from this file or its raw GitHub URL) |
| `zenodo-harvest-spm.gxwf.yml` | Source ([gxformat2](https://github.com/galaxyproject/gxformat2)); edit this one |
| `zenodo-harvest-spm-tests.yml` | `planemo test` cases for the `.ga` |

After editing the source, regenerate the `.ga` and test it:

```bash
gxwf-to-native zenodo-harvest-spm.gxwf.yml zenodo-harvest-spm.ga
planemo test --docker --conda_auto_install --extra_tools ../../tools/zenodo_nomad_harvest zenodo-harvest-spm.ga
```

The `.ga` uses the tools' short IDs. Once the tools are in the Tool Shed, re-export it from the Galaxy
server so it has their full Tool Shed IDs.

Author: Rubel Mozumder ([FAIRmat](https://www.fairmat-nfdi.eu/fairmat/))
