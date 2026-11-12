```yaml
Instrument:
  hardware:
    model: null
    model/@version: null
    name: Nanonis
    vendor: Nanonis
  software:
    model: null
    model/@version: "2"
    name: Nanonis
    vendor: Nanonis
Sample:
  History:
    Note:
      description: "The crystals were cleaved at approximately 100 K under ultra-high vacuum (UHV) conditions (P<10−10 mbar) just before being transferred to a low-temperature scanning tunneling microscope."
  Sample_component:
    chemical_formula: Bi2Pd
    description: "Bi-terminated surface of β-Bi2 Pd(001) crystals. The sample is superconducting below a critical temperature of 4.5 K, exhibiting a superconducting gap of ∆sample = 0.775 meV at this temperature."
    name: "β-Bi2 Pd single crystal"
  name: "β-Bi2 Pd single crystal"
default: current_forward
definition: NXstm
experiment_description: "Scanning tunneling microscopy was used to investigate the quasiparticle interference (QPI) and Yu-Shiba-Rusinov (YSR) states on the Bi-terminated surface of β-Bi2 Pd, revealing anisotropic modulations attributed to interband Cooper pair formation."
experiment_technique: STM
identifier_collection: https://doi.org/10.5281/zenodo.17202182
identifier_experiment: T190613_004.sxm
scan_mode: constant current
User:
  email:
    - operator@example.org
  name:
    - Example Operator
citeID:
  author:
    - Pascual, Jose Ignacio
  url:
    - https://zenodo.org/records/17202182
  doi:
    - 10.5281/zenodo.17202182
  description:
    - 'OPEN DATASET of the publication "Revealing band-hybrid Cooper pairs on the surface of a superconductor with spin-orbit coupling" including Scanning Tunneling Microscopy and Spectroscopy results.'
```

```yaml
evidence:
  Sample.History.Note.description:
    file: QPI_visualization_of_Interband_Pairing-SM.txt
    quote: "The crystals were cleaved at approximately 100 K under ultra-high vacuum (UHV) conditions (P<10−10 mbar) just before being transferred to a low-temperature scanning tunneling microscope"
  Sample.Sample_component.chemical_formula:
    file: QPI_visualization_of_Interband_Pairing-SM.txt
    quote: "The β−Bi2 Pd substrate becomes superconducting below a critical temperature of 4.5 K"
  Sample.Sample_component.description:
    file: QPI_visualization_of_Interband_Pairing-SM.txt
    quote: "The β−Bi2 Pd substrate becomes superconducting below a critical temperature of 4.5 K, exhibiting a superconducting gap of ∆sample = 0.775 meV at this temperature."
  Sample.Sample_component.name:
    file: QPI_visualization_of_Interband_Pairing-SM.txt
    quote: "The experiments were conducted on the clean Bi-terminated surface of β-Bi2 Pd(001) samples"
  Sample.name:
    file: QPI_visualization_of_Interband_Pairing-SM.txt
    quote: "The experiments were conducted on the clean Bi-terminated surface of β-Bi2 Pd(001) samples"
  experiment_description:
    file: QPI_visualization_of_Interband_Pairing-AcceptedMS.txt
    quote: "Here, we present evidence of the formation of Cooper pairs between spin-nondegenerate helical surface bands in the superconductor β-Bi2 Pd."
```
