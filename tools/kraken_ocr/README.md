# Kraken-OCR Suite

## Configuration for Galaxy Administrators

### Enabling of GPU Processing for Segmentation and OCR

By default, Kraken performs segmentation and OCR on the CPU. 
GPU execution requires setting the `KRAKEN_DEVICE` environment variable.

Allowed options:
- `cuda:0`, `cuda:1`,... -> use a specific CUDA device
- `cpu` -> use the CPU

Below is an example `job_conf.xml`:

```xml
<job_conf>
  <plugins>
    <plugin 
      id="local" 
      type="runner"
      load="galaxy.jobs.runners.local:LocalJobRunner"/>
  </plugins>

  <destinations default="local">
    <destination id="local" runner="local"/>

    <destination id="gpu" runner="local">
      <env id="KRAKEN_DEVICE">cuda:0</env>
    </destination>
  </destinations>

  <tools>
    <tool id="kraken_segment" destination="gpu"/>
    <tool id="kraken_ocr" destination="gpu"/>
  </tools>
</job_conf>
```

Note, GPU processing is only supported for kraken `segment` and `ocr`.

## Test-data Sources & Attribution

| File | Source | Modifications | License |
|---|---|---|---|
| `input.jpg` | From [`kraken/tests/resources/input.jpg`](https://github.com/mittagessen/kraken/blob/6ccf8d9a1c39d8bdc98d940e9764559c42344095/tests/resources/input.jpg) | None | Apache-2.0 |
| `input.tiff` | Derived from `input.jpg` | Converted to TIFF | Apache-2.0 |
| `input_2pages.pdf` | Derived from `input.jpg` and [`arabic.webp`](https://github.com/mittagessen/kraken/blob/37c3016f749d15715e8e17a8747d33f044fcb4b7/tests/resources/arabic.webp) | Grayscaled and combined into a two-page PDF | Apache-2.0 |
| `input_2pages.tiff` | Derived from `input_2pages.pdf` | Converted into a two-page TIFF | Apache-2.0 |
| `input_binarized.png` | Derived from `input.jpg` | Binarized with Kraken | Apache-2.0 |
| `input_binarized_seg_alto.xml` | Derived from `input_binarized.png` | Generated with Kraken 7.0.2 | Apache-2.0 |
| `input_binarized_seg_page.xml` | Derived from `input_binarized.png` | Generated with Kraken 7.0.2 | Apache-2.0 |
| `input_line.tiff` | Derived from [`OCR-kraken-models/example-snippest/sophoclesplaysa05campgoog_0177_19.png`](https://github.com/AjaxMultiCommentary/OCR-kraken-models/blob/57c18c88057fb3e16d17235500e608dac1c97566/example-snippets/sophoclesplaysa05campgoog_0177_19.png) | Converted to TIFF | CC-BY-4.0 |
| `lines_2pages.pdf` | Derived from `input_line.tiff` and [`OCR-kraken-models/example-snippets/sophokle1v3soph_0140_44.png`](https://github.com/AjaxMultiCommentary/OCR-kraken-models/blob/5c45ca7b1a66ac4115826bfaeb3f0476b2b1702a/example-snippets/sophokle1v3soph_0140_44.png) | Combined into a two-page PDF | CC-BY-4.0 |

### CC BY attribution for OCR-kraken-models example snippets

The CC-BY-4.0 source material was produced by the
[Ajax Multi-Commentary project](https://github.com/AjaxMultiCommentary/OCR-kraken-models/blob/e187951d0ebc669929c866f8a43bc303e849ff45/README.md).
Contributors according to the reposistory's README file are:
- Bruce Robertson (Mount Allison University)
- Sven Najem-Meyer (EPFL)
- Matteo Romanello (UNIL)