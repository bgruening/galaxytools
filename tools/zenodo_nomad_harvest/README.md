# zenodo_nomad_harvest

Galaxy tools that harvest scanning probe microscopy (SPM) data from [Zenodo](https://zenodo.org),
convert the raw files to [NeXus](https://www.nexusformat.org) and prepare
[NOMAD](https://nomad-lab.eu) projects. They are the steps of the
[Zenodo harvest workflow](../../workflows/zenodo_harvest/), which the `zenodo_harvest` Action of
[nomad-eosc-galaxy-actions](https://github.com/FAIRmat-NFDI/nomad-eosc-galaxy-actions) starts.

A file that cannot be handled is recorded with its reason; the job itself still succeeds.

## Tools, in workflow order

| Step | Tool | Input | Output |
|---|---|---|---|
| 1 | `zenodo_download_unzip` | `work_list.json` (from NOMAD); optional files already in Galaxy | `manifest.json`, raw files, one context bundle per raw file, context files |
| 2 | `nexus_eln_builder` | One context bundle | `eln_data.yaml`, ELN report (fields not found) |
| 3 | `pynxtools_dataconverter` | One raw file + its ELN | NeXus file (`.nxs`), status JSON, converter log |
| 4 | `nexus_report_builder` | Manifest and the outputs of steps 1–3 | `results.json`, one NOMAD project zip per record, logs per record, HTML summary |

- **Step 1** downloads the files of each record (MD5 checked) and extracts archives safely (no paths
  outside the folder, size limits). It gives each file a role: `raw`, `context` or `ignored`, from the
  work list's extension lists. Raw files the NOMAD side already knows (`known_sha256`) are skipped.
  The text of PDF context files is read with `pdftotext` into the bundle's `pdf_text`.
- **Step 2** fills the ELN from the record metadata and the text of the context files (readme, notes).
  It never guesses a value.
- **Step 3** runs `pynxtools` with the `spm` reader of
  [pynxtools-spm](https://github.com/FAIRmat-NFDI/pynxtools-spm) in a container (`macros.xml`).
- **Step 4** joins everything by file checksum (SHA-256) and lays out each project zip per Zenodo
  version: one folder per converted raw file and a `Complementaries` folder for the other files.

Steps 2 and 3 run once per raw file (Galaxy maps them over the collections). Raw files and ELNs share
element identifiers (`<record id>-<sha256 prefix>-<file name>`), so Galaxy pairs them.

## LLM ELN (optional workflow branch, on by default)

`nexus_llm_input` writes, per raw file, the text files of one
[LLM Hub](https://github.com/bgruening/galaxytools/tree/master/tools/llm_hub) job:

- `eln_instructions.txt`: fixed instructions per definition, from `eln_instructions/`;
- `base_eln.txt`: the rule ELN, with the fields the LLM may fill set to `null`;
- the context texts (PDF text included) and the `.sxm` header.

`eln_rules.yaml` holds the allowed `scan_mode` values per definition and the fields the LLM may
fill. The tool stops if an instructions file misses one of them, so the two stay in step.

`nexus_eln_llm_merge` checks the LLM answer and merges it into the rule ELN. A value is kept only
when its quote is in a source file and it passes the field's checks (allowed values, exact values,
no copied examples). Every kept, dropped or cleaned value is in its report. Without an answer the
rule ELN is the output. In the workflow, `llm_hub` runs between the two only when `use llm` is true.

## Exchange with NOMAD

- In: `work_list.json`: records, files to download, limits, extension lists, known raw files.
- Out: `manifest.json` (every archive and file, with role, action and reason), `results.json`
  (conversion status per raw file, by SHA-256) and the project zips.

Top-level keys (examples in `test-data/`):

| File | Keys |
|---|---|
| `work_list.json` | `schema_version`, `run_id`, `user`, `limits`, `raw_extensions`, `context_extensions`, `ignore_patterns`, `records` (each: record metadata, `known_sha256`, `download`) |
| `manifest.json` | `schema_version`, `run_id`, `tool_version`, `archives`, `files`, `download_failures` |

## Tests

```bash
planemo test --docker --conda_auto_install tools/zenodo_nomad_harvest
```

Which test files each tool uses, and what they check: [`test-data/README.md`](test-data/README.md).
