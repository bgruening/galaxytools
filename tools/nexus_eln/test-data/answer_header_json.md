Hand-written answer for the header JSON test: header quotes for header fields are kept; a key
alone and a sample name from the header are dropped.

```yaml
scan_mode: constant current
Instrument:
  hardware:
    model/@version: "2.1"
  software:
    model: Generic 4
Sample:
  name: 2H-NbSe2
```

```yaml
evidence:
  scan_mode:
    file: 7_T190613_004_header.json
    quote: '"Z-Controller::Controller status": "ON",'
  Instrument.software.model:
    file: 7_T190613_004_header.json
    quote: '"NanonisMain::SW Version": "Generic 4",'
  Instrument.hardware.model/@version:
    file: 7_T190613_004_header.json
    quote: '"Bias::Bias (V)":'
  Sample.name:
    file: 7_T190613_004_header.json
    quote: "2H-NbSe2"
```
