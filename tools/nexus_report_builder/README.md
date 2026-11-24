# nexus_report_builder

Joins the manifest of `zenodo_nomad_download` with the converter's outputs by file checksum (SHA-256).
Writes `results.json` (conversion status per raw file), one project zip per record for
[NOMAD](https://nomad-lab.eu), logs per record and an HTML summary. Each project zip has one folder per converted
raw file and a `Complementaries` folder for the other files.

## Tests

| Files | Checks |
|---|---|
| `manifest.json`, `status__*`, `nexus__*`, `convlog__*`, `eln_report__*`, `context__*`, `eln_data.yaml` | One converted and two failed raw files of one record; a file without status is `no_result` |
| `manifest_empty.json` | Empty run: empty results, no zips or logs |

File names with `__<record>-<hash>-<name>` are the element identifiers Galaxy gives in the workflow.
