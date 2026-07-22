# get example pretrained ocr models from various sources

```bash
ZENODO_URL=https://zenodo.org/record
FOLDER=files
MODEL_FILENAME=catmus-print-fondue-large.mlmodel
DOI_SUFFIX=10592716
wget -c -t 3 -O $MODEL_FILENAME $ZENODO_URL/$DOI_SUFFIX/$FOLDER/$MODEL_FILENAME

GITHUB_URL=https://raw.githubusercontent.com
REPO=AjaxMultiCommentary/OCR-kraken-models
GIT_HASH=e187951d0ebc669929c866f8a43bc303e849ff45
FOLDER=kraken-models/greek-english_porson_sophoclesplaysa05campgoog
MODEL_FILENAME=greek-english_porson_sophoclesplaysa05campgoog.mlmodel
wget -c -t 3 -O $MODEL_FILENAME $GITHUB_URL/$REPO/$GIT_HASH/$FOLDER/$MODEL_FILENAME
```
    