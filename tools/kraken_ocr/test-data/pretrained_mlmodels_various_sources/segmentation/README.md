# get example pretrained segmentation models from various sources

```bash
BASE_URL=https://raw.githubusercontent.com
REPO=michaelscho/bdd-segmentation-regions
GIT_HASH=b818488376bc635f1a121427979b33dfc11638df
MODEL_FILENAME=bdd-segmentation-regions-1.0.mlmodel
wget -c -t 3 -O $MODEL_FILENAME $BASE_URL/$REPO/$GIT_HASH/$MODEL_FILENAME
```
