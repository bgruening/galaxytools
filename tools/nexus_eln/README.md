# nexus_eln

Three tools that build the ELN (`eln_data.yaml`) the NeXus converter reads besides a raw file. They
share code and rules, so they are one repository.

| Tool | Input | Output |
|---|---|---|
| `nexus_eln_builder` | Context bundle; optional raw file | Rule ELN, and a report of the fields not found |
| `nexus_llm_input` | Bundle, rule ELN; optional raw file | The text files of one [LLM Hub](https://github.com/bgruening/galaxytools/tree/master/tools/llm_hub) job |
| `nexus_eln_llm_merge` | Bundle, rule ELN, answer; optional raw file | ELN with the checked values, and a report |

- `nexus_eln_builder` fills the ELN from the record metadata, the context files and the raw file
  header. It never guesses a value.
- `nexus_llm_input` writes the instructions for the definition (`eln_instructions/`), the rule ELN
  with the fields to fill set to `null`, the context texts and the raw file header.
- `nexus_eln_llm_merge` keeps an answer's value only when its quote is in a source file and it passes
  the field's checks. Without an answer the rule ELN is the output. The answer format is in its help.

`eln_rules.yaml` holds the fields to fill and the allowed values per definition. `nexus_llm_input`
stops if an instructions file misses one of them, so the two stay in step.

## Tests

| Tool | Files | Checks |
|---|---|---|
| builder | `bundle_with_context.json` | Values from the record and context files; the user comes from the work list |
| builder | `bundle_without_context.json` | No context, no user: missing fields are reported |
| builder | `bundle_with_context.json`, `const_dos_No14_003.sxm` | An STM `.sxm` header gives technique, scan mode and software |
| builder | `bundle_afm.json`, `afm_header_only.sxm` | An AFM `.sxm` header gives AFM; the technique parameter still wins |
| builder | `bundle_with_context.json`, `broken.sxm` | An unreadable header falls back to the extension |
| LLM input | `bundle_17202182.json`, `eln_17202182.yaml`, `T190613_004_header_only.sxm` | STM record [17202182](https://zenodo.org/records/17202182) (CC-BY-4.0): header, readme and PDF text |
| LLM input | `bundle_afm.json`, `eln_afm.yaml`, `afm_header_only.sxm` | AFM: readme and header; `scan_mode` left to the LLM |
| LLM input | `bundle_without_context.json`, `eln_sts.yaml` | STS without context or raw file |
| LLM input | `bundle_without_context.json`, `eln_no_technique.yaml` | No definition: only the base ELN |
| merge | `answer_afm_checks.md` (hand-written) | Each field hits one check; quotes with typographic spaces and hyphens or a missing space are still found |
| merge | `answer_mistral_medium_3_5.md` | Real answer of Mistral Medium 3.5 (LLM Hub, usegalaxy.eu): names made ASCII; a header key alone is no evidence |
| merge | no answer, `answer_skipped.json`, `answer_no_yaml.md` | LLM off, skipped step, unreadable answer: the rule ELN is the output |

`eln_17202182.yaml` and `eln_afm.yaml` are rule ELNs: they have no `Sample` part, because the builder
never guesses it.
