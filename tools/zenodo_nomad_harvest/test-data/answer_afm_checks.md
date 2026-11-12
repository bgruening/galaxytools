Hand-written answer for the checks of nexus_eln_llm_merge; each field hits one rule.
Some quotes use a no-break space, a no-break hyphen or a missing space, as real models write them.

```yaml
experiment_technique: AFM
definition: NXafm
experiment_description: NaCl islands were imaged in tapping mode.
scan_mode: Tapping Mode
Instrument:
  hardware:
    model: Cypher ES
    model/@version: "3.3"
  software:
    model: Nanonis 5.1.0
Sample:
  name: 2H-MoS2(001) flake
  History:
    Note:
      description: NaCl evaporated at 300 K onto Cu(111).
  Sample_component:
    name: NaCl
    chemical_formula: NaCl island
    description: NaCl islands on Cu(111).
User:
  name:
    - Invented Person
```

```yaml
evidence:
  experiment_description:
    file: eln_instructions.txt
    quote: "Never guess or use outside knowledge."
  scan_mode:
    file: 5_readme.txt
    quote: "Measured in tapping mode with a CypherES microscope"
  Instrument.hardware.model:
    file: 5_readme.txt
    quote: "with a Cypher ES microscope at room‑temperature"
  Instrument.hardware.model/@version:
    file: 5_readme.txt
    quote: "Measured in tapping mode with a Cypher ES microscope"
  Sample.name:
    file: 5_readme.txt
    quote: "AFM data of NaCl islands on Cu(111)."
  Sample.History.Note.description:
    file: 5_readme.txt
    quote: "The NaCl was evaporated at 300 K onto Cu(111)."
  Sample.Sample_component.name:
    file: notes.txt
    quote: "AFM data of NaCl islands on Cu(111)."
  Sample.Sample_component.chemical_formula:
    file: 5_readme.txt
    quote: "AFM data of NaCl islands on Cu(111)."
```
