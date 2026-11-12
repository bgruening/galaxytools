```yaml
experiment_technique: STM
definition: NXstm
identifier_collection: https://doi.org/10.5281/zenodo.17202182
identifier_experiment: T190613_004.sxm
experiment_description: A topographic image and lock‑in dI/dV map at +50 mV were recorded on the clean Bi‑terminated β‑Bi2Pd(001) surface to obtain quasiparticle‑interference patterns in the normal state, providing information on the band structure above the superconducting gap.
citeID:
  url:
  - https://zenodo.org/records/17202182
  author:
  - Pascual, Jose Ignacio
  doi:
  - 10.5281/zenodo.17202182
  description:
  - OPEN DATASET of the publication "Revealing band-hybrid Cooper pairs on the surface
    of a superconductor with spin-orbit coupling" including  Scanning Tunneling Microscopy
    and Spectroscopy results.
scan_mode: constant current
Instrument:
  hardware:
    model: JT-SPECS
    model/@version: null
  software:
    model: Nanonis
Sample:
  History:
    Note:
      description: The crystals were cleaved at approximately 100 K under ultra‑high vacuum (UHV) conditions (P<10−10 mbar) just before being transferred
  Sample_component:
    chemical_formula: Bi2Pd
    description: clean Bi-terminated surface of β‑Bi2Pd(001)
    name: beta-Bi2Pd(001) cleaved
  name: beta-Bi2Pd(001) cleaved
```

```yaml
evidence:
  experiment_description:
    file: 3_Information.txt
    quote: "Topography and lock-in dI/dV at 50 mV"
  Instrument.software.model:
    file: 6_T190613_004_header.txt
    quote: ":NANONIS_VERSION: 2"
  Instrument.hardware.model:
    file: 5_QPI_visualization_of_Interband_Pairing-SM.txt
    quote: "low-temperature scanning tunneling microscope (JT‑SPECS GmbH)"
  Sample.History.Note.description:
    file: 5_QPI_visualization_of_Interband_Pairing-SM.txt
    quote: "The crystals were cleaved at approximately 100 K under ultra‑high vacuum (UHV) conditions (P<10−10 mbar) just before being transferred"
  Sample.Sample_component.chemical_formula:
    file: 5_QPI_visualization_of_Interband_Pairing-SM.txt
    quote: "β‑Bi2Pd"
  Sample.Sample_component.description:
    file: 5_QPI_visualization_of_Interband_Pairing-SM.txt
    quote: "clean Bi-terminated surface of β‑Bi2Pd(001)"
  Sample.Sample_component.name:
    file: 5_QPI_visualization_of_Interband_Pairing-SM.txt
    quote: "β‑Bi2Pd"
  Sample.name:
    file: 5_QPI_visualization_of_Interband_Pairing-SM.txt
    quote: "β‑Bi2Pd"
```