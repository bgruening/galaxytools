# Reference measurement test data

The six `reference_speaker{1,2}_{param,logparam,points}.csv` fixtures were
generated locally with `quay.io/bgruening/new-fave:1.2.1` from 
`KY25A_1.mp3` audio and its corresponding transcript: `KY25A_1.TextGrid`.

```sh
docker run --rm \
    --user "$(id -u):$(id -g)" \
    --volume "$PWD:/data" \
    --workdir /data \
    quay.io/bgruening/new-fave:1.2.1 \
    fave-extract audio-textgrid KY25A_1.mp3 KY25A_1.TextGrid --speakers all
```

Each output CSV was split by `speaker_num`:

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
