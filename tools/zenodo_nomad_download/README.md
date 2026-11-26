# zenodo_nomad_download

Downloads the files of the records in a NOMAD harvest work list from [Zenodo](https://zenodo.org) (MD5 checked)
and extracts archives safely: no paths outside the folder, size limits. Each file gets a role, `raw`,
`context` or `ignored`, from the work list's extension lists. Raw files the caller already knows
(`known_sha256`) are skipped. The text of PDF context files is read with `pdftotext`.

A file that cannot be handled is recorded in the manifest with its reason; the job still succeeds.

## Input and outputs

- In: `work_list.json`, with the keys `schema_version`, `run_id`, `user`, `limits`, `raw_extensions`,
  `context_extensions`, `ignore_patterns` and `records` (each: record metadata, `known_sha256`,
  `download`). Optional: files already in Galaxy.
- Out: `manifest.json` (`schema_version`, `run_id`, `tool_version`, `archives`, `files`,
  `download_failures`), the raw files, one context bundle (JSON) per raw file, and the context files.
  Also fixed files the workflow uses in place of a failed job: `no_answer` (LLM Hub), `no_header`
  (header) and `no_eln` (ELN builder).
  Raw files and bundles share element identifiers (`<record id>-<sha256 prefix>-<file name>`).

## Tests

No network: the files of each work list are given as "prefetched" inputs.

| Files | Checks |
|---|---|
| `work_list.json`, `stm_archive.zip`, `traversal.zip`, `loose_raw.sxm`, `loose_readme.rtf`, `loose_mismatch.txt` | Raw, context, junk and known files; unsafe paths (`traversal.zip`), too large and MD5-mismatching files |
| `work_list_uncompressed_limit.json` | A zip above the uncompressed limit is not extracted |
| `work_list_empty.json` | No records: empty manifest and collections |
| `work_list_loose.json`, `context_unknown.csv` | Context files keep a fixed type, also with binary bytes |
| `work_list_pdf.json`, `context_sm_page2.pdf`, `context_broken.pdf`, `context_no_text.pdf` | PDF text in the bundle; a broken PDF and one without text are only reported |

`stm_archive.zip` and `traversal.zip` are ignored by the repository's `.gitignore` (`*.zip`) and were
added with `git add -f`. If one changes, update its `md5` and `size` in `work_list*.json`.
