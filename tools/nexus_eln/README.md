# nexus_eln

Three tools that build the ELN (`eln_data.yaml`) the NeXus converter reads besides a raw file. They
share code, so they are one repository.

The tools hold no domain knowledge. It comes from the pynxtools plugin as input files, e.g. the
`galaxy/` folder of [pynxtools-spm](https://github.com/FAIRmat-NFDI/pynxtools-spm):

- **rules** (`eln_rules.yaml`): definitions with their technique, raw file extensions, header markers
  and allowed values, the fields to fill and their checks;
- **instructions** (`eln_instructions.md`): one section `# <name>` per LLM job part: `Detection`,
  `Common`, one per NeXus definition (`NXstm`), `Output`;
- **header**: the raw file's metadata as JSON, from the plugin's raw file header tool.

| Tool | Input | Output |
|---|---|---|
| `nexus_llm_input`, stage detect | Bundle, rules, instructions, header | The text files of the detection [LLM Hub](https://github.com/bgruening/galaxytools/tree/master/tools/llm_hub) job |
| `nexus_eln_builder` | Bundle, rules, header, detection answer (optional) | Rule ELN, and a report of the fields not found |
| `nexus_llm_input`, stage eln | Bundle, rule ELN, rules, instructions, header | The text files of the ELN job |
| `nexus_eln_llm_merge` | Bundle, rule ELN, answer, rules, instructions, header | ELN with the checked values, and a report |

- `nexus_eln_builder` fills the ELN from the record metadata and the context files. The technique
  comes from the first of: the tool parameter; the detection answer, if its quote is in a source
  file; the header markers (only one technique matches); the extension. It never guesses a value.
- `nexus_llm_input` writes only the job's sections of the instructions, the header, `record.txt`
  (record title and description) and the context texts; in the ELN job also the rule ELN with the
  fields to fill set to `null`.
- `nexus_eln_llm_merge` keeps an answer's value only when its quote is in a source file and it passes
  the field's checks. Without an answer the rule ELN is the output. The answer format is in its help.

`nexus_llm_input` stops if a section is missing or misses a technique, field or allowed value of the
rules, so the two files stay in step.

## Tests

`eln_rules.yaml` and `eln_instructions.md` in `test-data` are copies of pynxtools-spm's
`src/pynxtools_spm/galaxy/` files; keep them the same. The `*_header.json` files are outputs of
the `spm_raw_header` tool (`tools/pynxtools_spm`): `const_dos_header.json` and `flt_header.json`
(a pynxtools-spm Bruker `.FLT` test file) read with gwyddionpy; `T190613_004_header.json` and
`afm_header.json` the leading text of header-only `.sxm` files. `answer_detect_*.md` are
hand-written detection answers.

| Tool | Files | Checks |
|---|---|---|
| builder | `bundle_with_context.json` | Values from the record and context files; the user comes from the work list |
| builder | `bundle_without_context.json` | No context, no user: missing fields are reported |
| builder | `bundle_with_context.json` | The technique parameter wins over the extension; an unknown one fails |
| builder | `afm_header.json`, `flt_header.json` | A header marker (text and JSON header) wins over the extension; the conflict is reported |
| builder | `answer_detect_afm.md`, `answer_detect_sts.md` | A found quote wins; a scan mode of another technique is dropped |
| builder | `answer_detect_no_quote.md` | A quote in no source is dropped; the extension decides |
| LLM input | `bundle_17202182.json`, `eln_17202182.yaml`, `T190613_004_header.json` | STM record [17202182](https://zenodo.org/records/17202182) (CC-BY-4.0): header, readme and PDF text |
| LLM input | `bundle_afm.json`, `eln_afm.yaml`, `afm_header.json` | AFM: readme and header; `scan_mode` left to the LLM |
| LLM input | `bundle_without_context.json`, `eln_sts.yaml`, `no_header.json` | STS without context or header |
| LLM input | `bundle_without_context.json`, `eln_no_technique.yaml` | No definition: only the base ELN |
| LLM input | `const_dos_header.json` | Detection job: the `Detection` section, a gwyddionpy header (JSON), `record.txt` |
| LLM input | `eln_instructions_incomplete.md` | A missing section, or one that misses rule values, fails |
| merge | `answer_afm_checks.md` (hand-written) | Each field hits one check; quotes with typographic spaces and hyphens or a missing space are still found |
| merge | `answer_mistral_medium_3_5.md` | Real answer of Mistral Medium 3.5 (LLM Hub, usegalaxy.eu): names made ASCII; a header key alone is no evidence |
| merge | no answer, `answer_skipped.json`, `answer_no_yaml.md` | LLM off, skipped step, unreadable answer: the rule ELN is the output |
| merge | `answer_header_json.md` (hand-written), `const_dos_header.json` | Header JSON quotes fill header fields only; a JSON key alone is no evidence |

`eln_17202182.yaml` and `eln_afm.yaml` are rule ELNs: they have no `Sample` part, because the builder
never guesses it.
