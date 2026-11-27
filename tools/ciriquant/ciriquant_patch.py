#!/usr/bin/env python
"""Make the CIRIquant 1.1.3 bioconda package runnable under Python 3.

The bioconda build only received a partial 2to3 conversion, so the package
still contains Python-2-only constructs that break the main pipeline on every
Python 3 release:

``CIRIquant/utils.py``
    ``from distutils.version import LooseVersion`` -- distutils was removed in
    Python 3.12, so the import aborts the run during YAML config parsing.

``CIRIquant/circ.py``
    * ``sorted(..., cmp=by_circ, ...)`` -- the ``cmp`` keyword was removed in
      Python 3, raising TypeError while the output GTF is written.
    * ``int / int`` used where integers are required (``file.read()``,
      ``pysam get_overlap()``, ``[iter()] * n``) -> TypeError.
    * ``start / 500`` used to build and query the annotation bin index ->
      float bin keys never match the integer lookups, so every circRNA loses
      its host gene annotation.

This script copies the installed package plus its sibling ``libs`` directory
(holding ``CIRI2.pl``) into the job directory, patches the copies and leaves
the copy first on ``PYTHONPATH``.  Nothing in the conda environment is
modified.
"""

import argparse
import os
import shutil
import sys

# The samtools version check calls the Python-2 builtin cmp() on
# distutils.LooseVersion objects (a mixed str/int comparison). Every samtools
# shipped by bioconda satisfies the required 1.9 minimum, so the check is
# replaced by a no-op instead of being ported.
UTILS_REPLACEMENTS = [
    (
        "    from distutils.version import LooseVersion\n",
        "",
    ),
    (
        "    if version and cmp(LooseVersion(version), "
        "LooseVersion('1.9')) < 0:\n"
        "        raise ConfigError('samtools version too low, "
        "1.9 required')\n",
        "    # samtools >= 1.9 is guaranteed by the conda environment\n",
    ),
]

CIRC_REPLACEMENTS = [
    # load_fai(): byte offset must stay an integer (passed to file.read())
    (
        "shift_length = int(length) * int(line_length) / int(eff_length)",
        "shift_length = int(length) * int(line_length) // int(eff_length)",
    ),
    # denovo_worker(): pysam get_overlap() needs integer coordinates
    (
        "circ_id, junc_site = d['SN'], int(d['LN']) / 2",
        "circ_id, junc_site = d['SN'], int(d['LN']) // 2",
    ),
    # proc_denovo_bam(): [iter(x)] * n needs an integer chunk size
    (
        "chunk_size = max(500, len(header) / thread + 1)",
        "chunk_size = max(500, len(header) // thread + 1)",
    ),
    # index_annotation() / circRNA_attr(): 500 bp bins must be integer keys
    (
        "start_div, end_div = parser.start / 500, parser.end / 500",
        "start_div, end_div = parser.start // 500, parser.end // 500",
    ),
    (
        "start_div, end_div = circ.start / 500, circ.end / 500",
        "start_div, end_div = circ.start // 500, circ.end // 500",
    ),
]

# sorted(..., cmp=by_circ, key=...): the cmp keyword was removed in Python 3.
# The original line sorts the circRNAs of a chromosome by the by_circ
# ordering of the parser objects, so it is rewritten with functools.cmp_to_key.
CIRC_SORT_LINE = (
    "for circ_id in sorted(list(circ_info[chrom].keys()), "
    "cmp=by_circ, key=lambda x:circ_info[chrom][x]):"
)
CIRC_SORT_REPLACEMENT = (
    "for circ_id in sorted(list(circ_info[chrom].keys()), "
    "key=functools.cmp_to_key(lambda x, y: by_circ("
    "circ_info[chrom][x], circ_info[chrom][y]))):"
)


def apply_replacements(src, replacements, label):
    for old, new in replacements:
        if old not in src:
            raise RuntimeError("%s: could not find %r" % (label, old[:60]))
        src = src.replace(old, new)
    return src


def main():
    parser = argparse.ArgumentParser(
        description="Patch the installed CIRIquant package for Python 3."
    )
    parser.add_argument("--dest", default="patched_pkg")
    args = parser.parse_args()

    import CIRIquant

    pkg_dir = os.path.dirname(os.path.abspath(CIRIquant.__file__))
    parent = os.path.dirname(pkg_dir)
    pkg_name = os.path.basename(pkg_dir)

    dest = os.path.abspath(args.dest)
    if os.path.isdir(dest):
        # Already patched by a previous invocation in this job directory.
        print("patched package already present: %s" % dest)
        return 0

    os.makedirs(dest)
    shutil.copytree(pkg_dir, os.path.join(dest, pkg_name))
    # main.py locates CIRI2.pl through <parent-of-package>/libs/CIRI2.pl
    libs = os.path.join(parent, "libs")
    if os.path.isdir(libs):
        shutil.copytree(libs, os.path.join(dest, "libs"))

    patched_pkg = os.path.join(dest, pkg_name)

    utils_path = os.path.join(patched_pkg, "utils.py")
    with open(utils_path) as fh:
        utils_src = fh.read()
    with open(utils_path, "w") as fh:
        fh.write(apply_replacements(utils_src, UTILS_REPLACEMENTS, "utils.py"))

    circ_path = os.path.join(patched_pkg, "circ.py")
    with open(circ_path) as fh:
        circ_src = fh.read()

    if CIRC_SORT_LINE not in circ_src:
        raise RuntimeError("circ.py: could not find the sorted(cmp=...) line")
    circ_src = circ_src.replace(CIRC_SORT_LINE, CIRC_SORT_REPLACEMENT)
    circ_src = apply_replacements(circ_src, CIRC_REPLACEMENTS, "circ.py")

    # functools is needed for the cmp_to_key based sort rewrite above
    circ_src = circ_src.replace(
        "import re\n",
        "import re\nimport functools\n",
        1,
    )

    with open(circ_path, "w") as fh:
        fh.write(circ_src)

    print("patched CIRIquant written to %s" % dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
