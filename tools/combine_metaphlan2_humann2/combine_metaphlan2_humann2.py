#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse


def extract_clade_abundance(metaphlan2_fp):
    clade_abund = {}
    with open(metaphlan2_fp, "r") as metaphlan2_f:
        for line in metaphlan2_f.readlines():
            if line.find("g__") == -1:
                continue

            split_line = line[:-1].split("\t")
            taxo = split_line[0]
            abundance = split_line[1]

            genus = taxo[(taxo.find("g__") + 3) :]
            if genus.find("|") != -1:
                genus = genus[: (genus.find("|"))]
            clade_abund.setdefault(genus, {"abundance": 0, "species": {}})
            if taxo.find("t__") != -1:
                continue
            elif taxo.find("s__") != -1:
                species = taxo[(taxo.find("s__") + 3) :]
                clade_abund[genus]["species"].setdefault(species, abundance)
            else:
                clade_abund[genus]["abundance"] = abundance
    return clade_abund


def compute_overall_abundance(humann2_fp):
    overall_abundance = 0
    with open(humann2_fp, "r") as humann2_f:
        for line in humann2_f.readlines():
            if line.find("|") != -1 or line.startswith("#"):
                continue
            split_line = line[:-1].split("\t")
            overall_abundance += float(split_line[1])
    return overall_abundance


def format_characteristic_name(name):
    formatted_n = name
    formatted_n = formatted_n.replace("/", " ")
    formatted_n = formatted_n.replace("-", " ")
    formatted_n = formatted_n.replace("'", "")
    if formatted_n.find("(") != -1 and formatted_n.find(")") != -1:
        open_bracket = formatted_n.find("(")
        close_bracket = formatted_n.find(")") + 1
        formatted_n = formatted_n[:open_bracket] + formatted_n[close_bracket:]
    return formatted_n


def combine_metaphlan2_humann2(args):
    clade_abund = extract_clade_abundance(args.metaphlan2_fp)
    overall_abund = compute_overall_abundance(args.humann2_fp)

    with open(args.output_fp, "w") as output_f:
        s = "genus\tgenus_abundance\tspecies\tspecies_abundance\t"
        s = "%s\t%s_id\t%s_name\t%s_abundance\n" % (s, args.type, args.type, args.type)
        output_f.write(s)
        with open(args.humann2_fp, "r") as humann2_f:
            for line in humann2_f.readlines():
                if line.find("|") == -1:
                    continue

                split_line = line[:-1].split("\t")
                abundance = 100 * float(split_line[1]) / overall_abund
                annotation = split_line[0].split("|")
                charact = annotation[0].split(":")
                charact_id = charact[0]
                char_name = ""
                if len(charact) > 1:
                    char_name = format_characteristic_name(charact[-1])
                taxo = annotation[1].split(".")

                if taxo[0] == "unclassified":
                    continue
                genus = taxo[0][3:]
                species = taxo[1][3:]

                if genus not in clade_abund:
                    print("no %s found in %s" % (genus, args.metaphlan2_fp))
                    continue
                if species not in clade_abund[genus]["species"]:
                    print(
                        "no %s found in %s for % s"
                        % (species, args.metaphlan2_fp, genus)
                    )
                    continue

                s = "%s\t%s\t" % (genus, clade_abund[genus]["abundance"])
                s += "%s\t%s\t" % (species, clade_abund[genus]["species"][species])
                s += "%s\t%s\t%s\n" % (charact_id, char_name, abundance)
                output_f.write(s)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--humann2_fp", required=True)
    parser.add_argument("--metaphlan2_fp", required=True)
    parser.add_argument("--output_fp", required=True)
    parser.add_argument("--type", required=True, choices=["gene_families", "pathways"])
    args = parser.parse_args()

    combine_metaphlan2_humann2(args)
