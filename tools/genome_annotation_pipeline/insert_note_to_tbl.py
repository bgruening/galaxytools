#!/usr/bin/env python

"""
    Insert note tag into the ncbi tbl file
"""

import argparse
import copy
import os
import pprint
import shutil
import sys
import tempfile

tbl_file = "STAFG.tbl"
discrep_file = "discrep.report"


def test():
    import filecmp

    tmp = tempfile.NamedTemporaryFile(delete=False)
    main("./tests/STAFG.tbl", "./tests/STAFG_discrep.report", tmp.name)
    assert filecmp.cmp("./tests/STAFG_result.tbl", tmp.name)
    os.remove(tmp.name)

    tmp = tempfile.NamedTemporaryFile(delete=False)
    main("./tests/STAFG.tbl", "./tests/1_discrep.report", tmp.name)
    assert filecmp.cmp("./tests/1_result.tbl", tmp.name)
    os.remove(tmp.name)

    print "All tests passed."


def main(tbl_file, discrep_file, outfile):
    discrep_overlapping_list = list()
    start = False
    for line in open(discrep_file):
        line = line.strip()
        if not line:
            start = False
            continue
        if line.startswith("DiscRep_SUB:OVERLAPPING_CDS"):
            start = True
            continue
        if start:
            ### typical line: STAFG:CDS	putative Precorrin-6Y C(5,15)-methyltransferase	lcl|Seq1234:1866-2915	STAFG_0415
            discrep_id, name, seq_pos, iid = line.split("\t")
            temp_dict = {
                "discrep_id": discrep_id,
                "name": name,
                "seq_pos": seq_pos,
                "iids": [iid],
            }
            if (
                discrep_overlapping_list
                and discrep_overlapping_list[-1]["name"] == name
            ):
                discrep_overlapping_list[-1]["iids"].append(iid)
            else:
                discrep_overlapping_list.append(temp_dict)

    tmp_out = False
    trash = list()
    for group in discrep_overlapping_list:

        if tmp_out:
            tmp_in = tmp_out.name
        else:
            # if we start that loop
            tmp_in = tbl_file

        tmp_out = tempfile.NamedTemporaryFile(delete=False)
        trash.append(tmp_out.name)

        found_product = False
        note_line = False
        for line in open(tmp_in):
            tmp_out.write(line)
            line = line.strip()
            if not line:
                continue

            # typical line: 			product	hypothetical protein
            if line.startswith("product") and line == "product\t%s" % group["name"]:
                found_product = True
                continue

            if found_product:
                # typical line: 			protein_id	gnl|PBUF|STAFG_1021
                for iid in group["iids"]:
                    if line == "protein_id\tgnl|PBUF|%s" % iid:
                        group_exclude = copy.deepcopy(group["iids"])
                        group_exclude.remove(iid)
                        if len(group_exclude) == 1:
                            note_line = (
                                "			note	overlaps CDS with %s\n" % group_exclude[0]
                            )
                        elif len(group_exclude) == 2:
                            note_line = "			note	overlaps CDS with %s and %s\n" % (
                                group_exclude[0],
                                group_exclude[1],
                            )
                        else:
                            note_line = "			note	overlaps CDS with %s and %s\n" % (
                                ", ".join(group_exclude[:-1]),
                                group_exclude[-1],
                            )
                        found_product = False

                found_product = False
            if note_line and line.startswith("note\t"):
                tmp_out.write(note_line)
                note_line = False
        tmp_out.close()
    shutil.copyfile(tmp_out.name, outfile)
    for filepath in trash:
        os.remove(filepath)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract differential methylated regions."
    )

    parser.add_argument("--tbl", required=True, help="NCBI table file.")

    parser.add_argument("--discrep", required=True, help="NCBI discrepancy file.")

    parser.add_argument("-o", "--outfile", required=True)

    parser.add_argument(
        "-t", "--test", action="store_true", default=False, help="Run tests."
    )

    options = parser.parse_args()
    main(options.tbl, options.discrep, options.outfile)

    if options.test:
        test()
