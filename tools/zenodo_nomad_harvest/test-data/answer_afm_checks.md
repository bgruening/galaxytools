Hand-written answer for the checks of nexus_eln_llm_merge; each field hits one rule.

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
  name: 2H-MoS2(0001) flake
  History:
    Note:
      description: NaCl evaporated at room temperature onto Cu(111).
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
    quote: "Measured in tapping mode with a Cypher ES microscope."
  Instrument.hardware.model:
    file: 5_readme.txt
    quote: "Measured in tapping mode with a Cypher ES microscope."
  Instrument.hardware.model/@version:
    file: 5_readme.txt
    quote: "Measured in tapping mode with a Cypher ES microscope."
  Sample.name:
    file: 5_readme.txt
    quote: "AFM data of NaCl islands on Cu(111)."
  Sample.History.Note.description:
    file: notes.txt
    quote: "The NaCl was evaporated at room temperature onto the clean Cu(111) crystal."
  Sample.Sample_component.name:
    file: 5_readme.txt
    quote: "AFM data of NaCl islands on Cu(111)."
  Sample.Sample_component.chemical_formula:
    file: 5_readme.txt
    quote: "AFM data of NaCl islands on Cu(111)."
```
