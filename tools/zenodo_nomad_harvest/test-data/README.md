# Test data

The tools run in this order in the workflow; each tool's tests use the files of its group.
Galaxy tool tests cannot chain tools, so the output of one step is saved here as the input of the next.

```
work list → zenodo_download_unzip → nexus_eln_builder → nexus_llm_input → (llm_hub) → nexus_eln_llm_merge
          → pynxtools_dataconverter → nexus_report_builder
```

Run all tool tests:

```bash
planemo test --docker --conda_auto_install tools/zenodo_nomad_harvest
```

## 1. `zenodo_download_unzip`: work list → manifest, raw files, context bundles

No network: the files of each work list are given as "prefetched" inputs.

| Files | Checks |
|---|---|
| `work_list.json`, `stm_archive.zip`, `traversal.zip`, `loose_raw.sxm`, `loose_readme.rtf`, `loose_mismatch.txt`, `const_dos_No14_003.sxm`, `broken.sxm` | Raw, context, junk and known files; unsafe paths (`traversal.zip`), too large and MD5-mismatching files |
| `work_list_uncompressed_limit.json` | A zip above the uncompressed limit is not extracted |
| `work_list_empty.json` | No records: empty manifest and collections |
| `work_list_loose.json`, `context_unknown.csv` | Context files keep a fixed type, also with binary bytes |
| `work_list_pdf.json`, `context_sm_page2.pdf`, `context_broken.pdf`, `context_no_text.pdf` | PDF text in the bundle; a broken PDF and one without text are only reported |

`stm_archive.zip` and `traversal.zip` are ignored by `.gitignore` (`*.zip`) and were added with
`git add -f`. If one changes, update its `md5` and `size` in `work_list*.json` and `manifest.json`.

## 2. `nexus_eln_builder`: bundle → rule ELN

| Files | Checks |
|---|---|
| `bundle_with_context.json` | Values from the record and context files; the user comes from the work list |
| `bundle_without_context.json` | No context, no user: missing fields are reported |
| `bundle_with_context.json`, `const_dos_No14_003.sxm` | An STM `.sxm` header gives technique, scan mode and software |
| `bundle_afm.json`, `afm_header_only.sxm` | An AFM `.sxm` header (recorded channels) gives AFM; the technique parameter still wins |
| `bundle_with_context.json`, `broken.sxm` | An unreadable header falls back to the extension |

## 3. `nexus_llm_input`: bundle + rule ELN → LLM Hub input files

| Files | Checks |
|---|---|
| `bundle_17202182.json`, `eln_17202182.yaml`, `T190613_004_header_only.sxm` | STM record [17202182](https://zenodo.org/records/17202182) (CC-BY-4.0): header, readme and PDF text; `base_eln.txt` sets the fields to fill to `null` |
| `bundle_afm.json`, `eln_afm.yaml`, `afm_header_only.sxm` | AFM: readme and header; `scan_mode` left to the LLM |
| `bundle_without_context.json`, `eln_sts.yaml` | STS without context or raw file |
| `bundle_without_context.json`, `eln_no_technique.yaml` | No definition: only the base ELN |

`eln_17202182.yaml` and `eln_afm.yaml` are rule ELNs: they have no `Sample` part, because the rule
builder never guesses it. Filling it is the LLM's job.

## 4. `nexus_eln_llm_merge`: LLM answer → checked, merged ELN

Uses the bundles and rule ELNs of group 3.

| Files | Checks |
|---|---|
| `answer_afm_checks.md` (hand-written) | Each field hits one check: quote not found, from the instructions, no evidence, value not in quote, copied example, invalid formula, allowed value; quotes with typographic spaces and hyphens or a missing space are still found |
| `answer_mistral_medium_3_5.md` | Real answer of Mistral Medium 3.5 (LLM Hub, usegalaxy.eu) for record 17202182: names made ASCII; a quote that is only a header key is dropped |
| no answer, `answer_skipped.json`, `answer_no_yaml.md` | LLM off, skipped step, unreadable answer: the rule ELN is the output |

## 5. `pynxtools_dataconverter`: raw file + ELN → NeXus

| Files | Checks |
|---|---|
| `const_dos_No14_003.sxm`, `eln_data.yaml` | A Nanonis STM image converts to NXstm |
| `broken.sxm`, `raw_without_extension` | Recorded failures; the job stays green |
| `eln_no_technique.yaml` | `auto` without a definition is a recorded failure |

## 6. `nexus_report_builder`: manifest + statuses → results, project zips, logs

| Files | Checks |
|---|---|
| `manifest.json`, `status__*`, `nexus__*`, `convlog__*`, `eln_report__*`, `context__*`, `eln_data.yaml` | One converted and two failed raw files of one record; a file without status is `no_result` |
| `manifest_empty.json` | Empty run: empty results, no zips or logs |

File names with `__<record>-<hash>-<name>` are the element identifiers Galaxy gives in the workflow.

## Workflow tests

`workflows/zenodo_harvest/zenodo-harvest-spm-tests.yml` runs the whole chain with the LLM off, using
`work_list.json`, `work_list_empty.json` and `work_list_loose.json` with their files from group 1.
