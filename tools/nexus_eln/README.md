# nexus_eln

Three tools that build the ELN (`eln_data.yaml`) the NeXus converter reads besides a raw file. They
share code, so they are one repository.

The tools hold no domain knowledge. It comes from the pynxtools plugin as input files, e.g. the
`galaxy/` folder of [pynxtools-spm](https://github.com/FAIRmat-NFDI/pynxtools-spm):

- **rules** (`eln_rules.yaml`): definitions with their technique, raw file extensions and allowed
  values, the fields to fill and their checks;
- **instructions**: one text file per NeXus definition (`NXstm.txt`), as a list collection;
- **header**: the raw file's metadata as JSON, from the plugin's raw file header tool.

| Tool | Input | Output |
|---|---|---|
| `nexus_eln_builder` | Context bundle, rules | Rule ELN, and a report of the fields not found |
| `nexus_llm_input` | Bundle, rule ELN, rules, instructions, header | The text files of one [LLM Hub](https://github.com/bgruening/galaxytools/tree/master/tools/llm_hub) job |
| `nexus_eln_llm_merge` | Bundle, rule ELN, answer, rules, instructions, header | ELN with the checked values, and a report |

- `nexus_eln_builder` fills the ELN from the record metadata and the context files; the technique
  comes from the raw file's extension, as listed in the rules. It never guesses a value.
- `nexus_llm_input` writes the instructions for the definition, the rule ELN with the fields to fill
  set to `null`, the context texts and the header.
- `nexus_eln_llm_merge` keeps an answer's value only when its quote is in a source file and it passes
  the field's checks. Without an answer the rule ELN is the output. The answer format is in its help.

`nexus_llm_input` stops if an instructions file misses a field or an allowed value of the rules, so
the two stay in step.

## Tests

`eln_rules.yaml` and `eln_instructions/` in `test-data` are copies of pynxtools-spm's
`src/pynxtools_spm/galaxy/` files; keep them the same. The `*_header.json` files are outputs of
the `raw_file_header` tool (`tools/pynxtools_spm`): `const_dos_header.json` read with gwyddionpy;
`T190613_004_header.json` and `afm_header.json` the leading text of header-only `.sxm` files.

| Tool | Files | Checks |
|---|---|---|
| builder | `bundle_with_context.json` | Values from the record and context files; the user comes from the work list |
| builder | `bundle_without_context.json` | No context, no user: missing fields are reported |
| builder | `bundle_with_context.json` | The technique parameter wins over the extension; an unknown one fails |
| LLM input | `bundle_17202182.json`, `eln_17202182.yaml`, `T190613_004_header.json` | STM record [17202182](https://zenodo.org/records/17202182) (CC-BY-4.0): header, readme and PDF text |
| LLM input | `bundle_afm.json`, `eln_afm.yaml`, `afm_header.json` | AFM: readme and header; `scan_mode` left to the LLM |
| LLM input | `bundle_without_context.json`, `eln_sts.yaml`, `no_header.json` | STS without context or header |
| LLM input | `bundle_without_context.json`, `eln_no_technique.yaml` | No definition: only the base ELN |
| LLM input | `const_dos_header.json` | A gwyddionpy header as JSON; instructions named without extension; instructions that miss rule values fail |
| merge | `answer_afm_checks.md` (hand-written) | Each field hits one check; quotes with typographic spaces and hyphens or a missing space are still found |
| merge | `answer_mistral_medium_3_5.md` | Real answer of Mistral Medium 3.5 (LLM Hub, usegalaxy.eu): names made ASCII; a header key alone is no evidence |
| merge | no answer, `answer_skipped.json`, `answer_no_yaml.md` | LLM off, skipped step, unreadable answer: the rule ELN is the output |
| merge | `answer_header_json.md` (hand-written), `const_dos_header.json` | Header JSON quotes fill header fields only; a JSON key alone is no evidence |

`eln_17202182.yaml` and `eln_afm.yaml` are rule ELNs: they have no `Sample` part, because the builder
never guesses it.
