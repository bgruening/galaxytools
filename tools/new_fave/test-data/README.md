
# Reference measurement test data

The six `reference_speaker{1,2}_{param,logparam,points}.csv` fixtures were
generated locally with [`quay.io/bgruening/new-fave:1.2.1`](https://quay.io/repository/bgruening/new-fave?tab=tags) from 
`KY25A_1.mp3` audio and its corresponding transcript: `KY25A_1.TextGrid`
as follows:

## 1. Extract formant measurements and recoded TextGrid from audio and transcript:

```sh
docker run --rm \
    --user "$(id -u):$(id -g)" \
    --volume "$PWD:/data" \
    --workdir /data \
    quay.io/bgruening/new-fave:1.2.1 \
    fave-extract audio-textgrid KY25A_1.mp3 KY25A_1.TextGrid --speakers all
```
## 2. Split the output CSV by speaker number into the corresponding reference fixture:

Each output CSV in the `fave_results` directory was split by `speaker_num` (1 or 2):

```sh
for kind in param logparam points; do
    for speaker in 1 2; do
        awk -F, -v speaker="$speaker" '
            NR == 1 {
                for (i = 1; i <= NF; i++)
                    if ($i == "speaker_num") column = i
                print
                next
            }
            $column == speaker
        ' "fave_results/KY25A_1_${kind}.csv" \
          > "reference_speaker${speaker}_${kind}.csv"
    done
done
```

# Alignment and overlap fixtures

`KY25A_1_fave.TextGrid` is derived from `KY25A_1.TextGrid` by reordering
its four tiers from words/phones/words/phones to phones/words/phones/words
(the classic FAVE layout). Tier names, interval boundaries, and labels are
unchanged; both grids use `KY25A_1.mp3`. This is a format conversion, not a
new alignment produced by fave-align.

The original recording and grid already contain speech from two speakers
at the same time. The paired overlap tests use this unchanged input with
`exclude_overlaps` disabled and enabled, with optimization disabled in both
runs. In new-fave 1.2.1, 65 vowel measurements become 49: seven KY25A vowels
and nine IVR vowels are excluded. Tests check the removed interval IDs
explicitly in both runs and constrain the retained IDs in the points CSV,
without comparing floating-point formant measurements. The FAVE layout
test expects the same 65 interval IDs as the unfiltered original grid.

# Audio format fixtures

`KY25A_1.wav` and `KY25A_1.flac` are conversions of `KY25A_1.mp3`,
with the original sample rate and channel count preserved. The defaults test
uses WAV, the speakers configuration test uses FLAC, and other tests retain
MP3 coverage. Both conversions use the original `KY25A_1.TextGrid`.

Regenerate from this directory with FFmpeg:

```sh
ffmpeg -i KY25A_1.mp3 -c:a pcm_s16le KY25A_1.wav
ffmpeg -i KY25A_1.wav -c:a flac KY25A_1.flac
```
