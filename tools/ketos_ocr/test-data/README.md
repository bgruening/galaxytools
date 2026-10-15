# Ketos OCR test data

The `ketos_train_recognition.arrow` fixture was compiled with Kraken 7.1 from
[`170025120000003,0074-lite.xml`](https://github.com/mittagessen/kraken/blob/main/tests/resources/170025120000003%2C0074-lite.xml)
and its referenced JPEG. This is the same compact PAGE fixture used by Kraken's
upstream recognition training smoke test.

`ketos_train_model.safetensors` contains the best weights produced by a
one-epoch CPU run over that Arrow dataset using the small upstream smoke-test
VGSL specification:

```text
[1,12,0,1 Cr3,3,8 S1(1x0)1,3]
```

This is a deliberately tiny [Kraken VGSL network](https://kraken.re/5.2/vgsl.html)
chosen to keep the smoke tests fast:

- `[1,12,0,1]` defines the input as `[batch, height, width, channels]`: one
  12-pixel-high, variable-width (`0`) grayscale channel.
- `Cr3,3,8` applies a 3 × 3 convolution with ReLU activation (`r`) and eight
  output channels.
- `S1(1x0)1,3` collapses the height into the channel dimension. It splits
  dimension 1 (height) into `1 × 0`, where `0` means infer the remaining
  factor, leaves the size-1 part in dimension 1, and moves the inferred part
  to dimension 3 (channels). After the convolution, this changes each
  width-wise feature vector from height 12 × 8 channels to height 1 × 96
  channels.

The specification intentionally omits an output block: during recognition
training Ketos appends the CTC output layer sized to the alphabet found in the
training data. This minimal convolution-and-reshape network is suitable for
testing the training workflow, not for producing an accurate OCR model.

`ketos_train_validation_1.arrow` and `ketos_train_validation_2.arrow` contain
the first and last two records, respectively, of `ketos_train_recognition.arrow`.
Their metadata record counts are updated to two. These fixtures exercise explicit
validation with one or multiple files, and are also reused as two distinct training
inputs with split or explicit validation. They overlap the training fixture because
these are workflow smoke tests, not model-quality evaluations.

`ketos_train_resume.ckpt` is a full Kraken 7.1 / PyTorch Lightning 2.6.1
checkpoint after one epoch (epoch 0, global step 3), including optimizer and
training-loop state. It uses the same compact VGSL specification above and four
copies of `compile_input_line.png` / `compile_input_line.gt.txt`. The resume test
runs to a total of two epochs and checks that an epoch-1 checkpoint is exported.

The fresh-run checkpoint test also loads this fixture through `--load`, stages it
with a `.ckpt` suffix, and checks that training restarts at epoch 0 rather than
restoring the training-loop state.

To regenerate it, stage those four pairs in a temporary working directory as
`ground_truth_0.png` / `ground_truth_0.gt.txt` through
`ground_truth_3.png` / `ground_truth_3.gt.txt`, then run:

```sh
ketos --device cpu --workers 0 --threads 1 --deterministic --seed 42 train \
    --output model --weights-format safetensors --format-type path \
    --arch vgsl --spec '[1,12,0,1 Cr3,3,8 S1(1x0)1,3]' \
    --quit fixed --epochs 1 --freq 1 --no-augment --partition 0.75 \
    ground_truth_0.png ground_truth_1.png ground_truth_2.png ground_truth_3.png
```

Copy the resulting `model/checkpoint_00-*.ckpt` to `ketos_train_resume.ckpt`.
The relative input names intentionally match the wrapper's staging names:
Kraken restores the checkpoint's data configuration when resuming, so absolute
paths or paths into `test-data` would not work in a Galaxy job directory.

The checkpoint input uses `ftype="zip"` (a binary subtype) because PyTorch
checkpoints are ZIP containers. This keeps Galaxy from decompressing the upload:
Kraken needs the container intact to memory-map and restore it.

`ketos_train_resume_missing_data.ckpt` is derived from the resume fixture above,
with only its saved data configuration changed to use a missing compiled Arrow
dataset. The negative resume test supplies a valid binary dataset but expects the
saved path to be used, the file-opening warning to be printed, and the job to fail
with no training data. To regenerate this fixture in the Kraken environment:

```python
import torch

checkpoint = torch.load("test-data/ketos_train_resume.ckpt",
                        map_location="cpu", weights_only=False)
config = checkpoint["datamodule_hyper_parameters"]["data_config"]
config.format_type = "binary"
config.training_data = ["missing_resume_data/training.arrow"]
config.evaluation_data = []
config.test_data = []
torch.save(checkpoint, "test-data/ketos_train_resume_missing_data.ckpt")
```

`ketos_train_binary_resume.ckpt` is a one-epoch checkpoint from the compact VGSL
model trained on both validation Arrow fixtures, also used as explicit validation
inputs. Its saved paths are `training_0.arrow`, `training_1.arrow`,
`validation_0.arrow`, and `validation_1.arrow`. The binary resume Galaxy test
recreates those names and checks that training continues into epoch 1.

To regenerate it, stage the two validation Arrow fixtures under both pairs of
names in a temporary directory, write `validation_0.arrow` and
`validation_1.arrow` on separate lines in `evaluation_manifest.txt`, and run:

```sh
ketos --device cpu --workers 0 --threads 1 train \
    --output model --weights-format safetensors --format-type binary \
    --arch vgsl --spec '[1,12,0,1 Cr3,3,8 S1(1x0)1,3]' --batch-size 1 \
    --quit fixed --epochs 1 --freq 1 --no-augment \
    --evaluation-data evaluation_manifest.txt training_0.arrow training_1.arrow
```

Copy `model/checkpoint_00-*.ckpt` to `test-data/ketos_train_binary_resume.ckpt`.
The separate `tests/test_binary_resume.py` integration test generates fresh
checkpoints, removes the original job and its `.dat` inputs, and resumes using
new input paths. It covers both split and explicit validation with multiple
datasets. Run it with Galaxy's Python dependencies and Cheetah3 installed:

```sh
KETOS_TEST_EXECUTABLE=/path/to/kraken-environment/bin/ketos \
    python -m unittest discover -s tests -v
```

The Kraken environment must also contain TensorBoard. Without
`KETOS_TEST_EXECUTABLE`, only the form and command-rendering checks run.

`ketos_train_codec.json` is the character-to-label mapping from
`ketos_train_model.safetensors`, with an extra `~` character assigned label 36.
Label 0 is reserved for CTC blank. The extra character makes this codec differ
from the inferred training alphabet and the loaded model codec. It exercises
explicit JSON codecs for new VGSL training and loading with `resize=new`.
