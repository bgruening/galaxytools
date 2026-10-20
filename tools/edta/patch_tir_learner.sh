#!/usr/bin/env bash

# Workaround for two upstream TIR-Learner 3.0 bugs that abort EDTA's TIR module.
#
# EDTA 2.3.0 depends on TIR-Learner 3.0, which was written against pandas 1.x.
# With pandas >= 2 the TIR module dies and takes the whole EDTA run down with it:
#
#   1. check_TIR_TSD.py indexes row Series positionally (``x[0]``). pandas 2
#      removed that fallback, so the step aborts with ``KeyError: 0``.
#      Upstream report: https://github.com/oushujun/EDTA/issues/631
#   2. get_fasta_sequence.py slices the genome with ``x["start"]:x["end"]``
#      where the bounds are strings, aborting with
#      ``TypeError: slice indices must be integers``.
#      Upstream report: https://github.com/oushujun/EDTA/issues/651
#
# Both are patched in place until upstream ships fixes. The patch is best-effort:
# if the installed layout differs, the files are left untouched and EDTA runs
# unmodified, so this can never be the reason a job fails.

set -u

# The biocontainer installs into /usr/local, a conda env into $CONDA_PREFIX.
roots=()
if [ -n "${CONDA_PREFIX:-}" ]; then
    roots+=("${CONDA_PREFIX}")
fi
roots+=("/usr/local")

for root in "${roots[@]}"; do
    for dir in "${root}"/share/TIR-Learner*/bin; do
        [ -d "${dir}" ] || continue

        tsd="${dir}/check_TIR_TSD.py"
        if [ -f "${tsd}" ] && grep -q 'x\[[0-9]\+\]' "${tsd}"; then
            sed -i 's/x\[\([0-9]\+\)\]/x.iloc[\1]/g' "${tsd}" \
                && echo "Patched positional pandas indexing in ${tsd}"
        fi

        fasta="${dir}/get_fasta_sequence.py"
        if [ -f "${fasta}" ] && grep -q 'genome_SeqRecord\.seq\[x\[' "${fasta}"; then
            sed -i 's/genome_SeqRecord\.seq\[x\["start"\]: *x\["end"\]\]/genome_SeqRecord.seq[int(x["start"]):int(x["end"])]/g' "${fasta}" \
                && echo "Patched string slice indices in ${fasta}"
        fi
    done
done

# Never let a failed patch abort the EDTA run.
exit 0
