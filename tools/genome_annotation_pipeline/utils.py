#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import re
import sys


def change_according_reviewer(accession, note_line):
    accession = accession.strip()

    """
    putative D-malate dehydrogenase [decarboxylating] [gnl|PBUF|STVIR_0046:1-352] [gnl|PBUF|STVIR_0046: raw, aa len= 352]
    Saccharopine dehydrogenase [NAD(+), L-lysine-forming]
    Delete the brakets ...
    """
    words = []
    brackets_open = False
    for word in accession.split():
        word = word.strip()
        if word[0] == "[" and word[-1] == "]":
            continue
        if word[0] == "[":
            brackets_open = True
            continue
        elif word[-1] == "]" and brackets_open:
            brackets_open = False
            continue
        elif brackets_open:
            continue
        else:
            words.append(word)
    accession = " ".join(words)

    if (
        accession.strip() == ""
        or accession.lower().find("uncharacterized protein") != -1
    ):
        return ""

    """
        DiscRep_SUB:SUSPECT_PRODUCT_NAMES::Organelles not appropriate in prokaryote
        DiscRep_SUB:SUSPECT_PRODUCT_NAMES::12 features contains 'Chloroplast'
    """
    if accession.lower().find("chloroplast") != -1:
        return ""

    accession = accession.replace("(Fragments)", "")

    if accession.lower().find("truncated") != -1:
        return ""

    accession = accession.replace(
        "55.5 kDa and 49.5 kDa sporulation proteins", "sporulation protein"
    )

    """
        BRA0749/BS1330_II0742
    """
    accession = re.sub("\w*\d\d\d\d\/.*\d\d\d\d_\w*", "", accession)  # noqa W605

    """
        Meiotically up-regulated gene 72 protein -> skip it
    """
    if accession.lower().find("meiotically up-regulated gene") != -1:
        return ""

    accession = accession.split("OS=")[0]
    accession = re.sub("^[p,P]otential ", "", accession)
    accession = re.sub("^[p,P]robable ", "", accession)
    accession = accession.replace("sulpho", "sulfo")
    accession = accession.replace("FMRF-amide neuropeptides", "FMRF-amide neuropeptide")
    accession = accession.replace(" (Fragment)", "")
    accession = accession.replace("Staphylococcal ", "")

    # delete numbers and ids -- do not include locus_tags or database identifer in product names)
    accession = re.sub("\s\w*\d{3}[\w\d.]*", "", accession)  # noqa W605

    # contains 'Homolog' (use -like protein instead)
    # special case: Beige protein homolog 1
    accession = re.sub("homolog[\s\d]*$", "like protein", accession, re.IGNORECASE)  # noqa W605
    # normal use case
    accession = re.sub("homolog", "like protein", accession, re.IGNORECASE)

    # ends with 'like' (add protein to the end xxx-like protein)
    accession = re.sub("-like$", "-like protein", accession, re.IGNORECASE)

    # some accessions longer than 100 characters (use concise biological name, move long descriptions to a note)
    """
    accession.replace('SWI/SNF-related matrix-associated actin-dependent regulator of chromatin subfamily A member 3-like 3', 'SWI/SNF-related matrix-associated actin-dependent regulator')
    accession.replace('Dihydrolipoyllysine-residue acetyltransferase component of pyruvate dehydrogenase complex, mitochondrial', 'Dihydrolipoyllysine-residue acetyltransferase')
    accession.replace('1-(5-phosphoribosyl)-5-[(5-phosphoribosylamino)methylideneamino] imidazole-4-carboxamide isomerase', 'imidazole-4-carboxamide isomerase')
    accession.replace('Phosphatidylinositol-3,4,5-trisphosphate 3-phosphatase and dual-specificity protein phosphatase PTEN', 'Phosphatidylinositol-3,4,5-trisphosphate 3-phosphatase')
    accession.replace('Dihydrolipoyllysine-residue succinyltransferase component of 2-oxoglutarate dehydrogenase complex, mitochondrial', 'Dihydrolipoyllysine-residue succinyltransferase')
    accession.replace('SWI/SNF-related matrix-associated actin-dependent regulator of chromatin subfamily A member 3-like 1', 'SWI/SNF-related matrix-associated actin-dependent')
    accession.replace('Polyprenol-phosphate-mannose-dependent alpha-(1-2)-phosphatidylinositol mannoside mannosyltransferase', 'Alpha-D-mannose-alpha-(1-2)-mannosyltransferase')
    accession.replace('Dihydrolipoyllysine-residue succinyltransferase component of 2-oxoglutarate dehydrogenase complex', '2-oxoglutarate dehydrogenase complex component E2')
    accession.replace('UDP-N-acetylglucosamine--N-acetylmuramyl-(pentapeptide) pyrophosphoryl-undecaprenol N-acetylglucosamine transferase', 'Undecaprenyl-PP-MurNAc-pentapeptide-UDPGlcNAc GlcNAc transferase')
    accession.replace('Polyprenol-phosphate-mannose-dependent alpha-(1-2)-phosphatidylinositol mannoside mannosyltransferase', 'Alpha-D-mannose-alpha-(1-2)-mannosyltransferase')
    """

    # DiscRep_SUB:DISC_PRODUCT_NAME_TYPO::341 features contains 'putative probable', Replace with 'putative'
    accession = accession.replace("putative probable", "putative")
    accession = accession.replace("probable", "")
    accession = accession.replace("Uncharacterized", "")
    accession = accession.replace("Phosphohexose mutases", "Phosphohexose mutase")

    # DiscRep_SUB:DISC_PRODUCT_NAME_TYPO::4 features contains 'diacyglycerol', Replace with 'diacylglycerol'
    accession = accession.replace("diacyglycerol", "diacylglycerol")

    accession = accession.replace("putative", "", 1).replace("Putative", "", 1).strip()
    accession = accession.replace(
        "peptidase S8/S53 subtilisin kexin sedolisin", "peptidase"
    )
    accession = accession.replace("Alcohol", "alcohol")
    accession = accession.replace("Pro", "pro")
    accession = accession.replace("MFS_1", "MFS-1")
    accession = accession.replace("Endo", "endo")
    accession = accession.replace("Transcriptional", "transcriptional")
    accession = accession.replace("Xylan", "xylan")
    accession = accession.replace("Mandelate", "mandelate")
    accession = accession.replace("Phytanoyl", "phytanoyl")
    accession = accession.replace("Silent", "silent")
    accession = accession.replace("delta'", "delta")
    accession = accession.replace("Relaxase", "relaxase")
    accession = accession.replace("XRE", "xre")
    accession = accession.replace("Holliday", "holliday")
    accession = accession.replace("Membrane", "membrane")
    accession = accession.replace("Methyltransferase", "methyltransferase")
    accession = accession.replace("Epoxide", "epoxide")
    accession = accession.replace(
        "EA31 gene protein, phage lambda", "EA31 protein, phage lambda"
    )
    accession = accession.replace("disulphide", "disulfide")
    accession = accession.replace("utilising", "utilizing")
    accession = accession.replace("hypothetical membrane protein", "membrane protein")
    accession = accession.replace(
        "involving differentiation", "differentiation protein"
    )
    accession = accession.replace(
        "bifunctional protein (methylenetetrahydrofolate dehydrogenase and methenyltetrahydrofolate cyclohydrolase)",
        "bifunctional protein methylenetetrahydrofolate dehydrogenase and methenyltetrahydrofolate cyclohydrolase",
    )
    accession = accession.replace(
        "N- superfamily bifunctional DNA primase/polymerase",
        "bifunctional DNA primase/polymerase",
    )

    accession = accession.replace("accessroy", "accessory")
    accession = accession.replace("Putative", "")
    accession = accession.replace("putative", "")
    accession = accession.replace("Puative", "")
    beta_lactamase = re.compile(
        "Beta-lactamase class C and other penicillin binding proteins".lower(),
        re.IGNORECASE,
    )
    accession = beta_lactamase.sub("beta-lactamase class C", accession)

    accession = accession.replace(
        "NADH dehydrogenase I, B/C/D subunits", "NADH dehydrogenase I"
    )
    accession = accession.replace(
        "putative FAD-dependent pyridine nucleotide-disulphide oxidoreductase",
        "FAD-dependent pyridine nucleotide-disulfide oxidoreductase",
    )
    accession = accession.replace("dyhydrogenase", "dehydrogenase")
    accession = accession.replace(
        "Type III effector Hrp-dependent outers",
        "Type III effector Hrp-dependent outer protein",
    )
    accession = accession.replace("oxidoreducatse", "oxidoreductase")
    accession = accession.replace("putaive", "")
    accession = accession.replace("puitative", "")
    accession = accession.replace(
        "iron-sulphur-binding reductase", "iron-sulfur-binding reductase"
    )
    accession = accession.replace(
        "ubiquinol-cytochrome c reductase iron-sulphur",
        "ubiquinol-cytochrome c reductase iron-sulfur",
    )
    accession = accession.replace("trancsriptional", "transcriptional")
    accession = accession.replace("pXO2-42/BXB0045/GBAA_pXO2_0045", "")
    accession = accession.replace("/y1452/YP_2654", "")
    accession = accession.replace("pXO2-48/BXB0055/GBAA_pXO2_0055", "")
    accession = accession.replace("Catabolite gene activator", "cAMP receptor protein")

    if not note_line:
        accession = accession.replace("pCQ3_82", "hypothetical protein")
        accession = accession.replace("pCQ3_23", "hypothetical protein")
        accession = accession.replace("pCQ3_24", "hypothetical protein")
        accession = accession.replace(
            "sp|Q08972|NEW1_YEAST", "[NU+] prion formation protein 1"
        )
        accession = accession.replace(
            "sp|Q06449|PIN3_YEAST", "[PSI+] inducibility protein 3"
        )
        accession = accession.replace(
            "sp|Q9P6P9|PDK_SCHPO", "Pyruvate dehydrogenase kinase"
        )
        accession = re.sub(
            "sp\|\w{6}\|\w+_\w+", "", accession  # noqa W605
        )  # replace names like sp|P51831|FABG_BACSU

        # filter out swissprot ID's sp|O69873|GLND_STRCO
        # if accession.find('sp|') != -1 and accession.count('|') >= 2:
        #    return ''

    accession = accession.replace(
        "bacterial capsule synthesis protein PGA_cap",
        "bacterial capsule synthesis protein PGA cap",
    )
    mandelate = re.compile(
        "mandelate racemase/muconate lactonizing enzyme:Mandelate racemase/muconate lactonizing enzyme".lower(),
        re.IGNORECASE,
    )
    accession = mandelate.sub("beta-lactamase class C", accession)

    accession = accession.replace(
        "conserved hypothetical protein", "hypothetical protein"
    )

    accession = accession.replace("LOW QUALITY PROTEIN: ", "")
    if accession in ["membrane protein SC4H802", "membrane protein SC1C221c"]:
        accession = "membrane protein"
    if accession == "GF18879":
        accession = "hypothetical protein"
    if accession == "LOW QUALITY PROTEIN: hypothetical protein SSLG_04259":
        accession = "hypothetical protein"
    if accession == "F420-0--gamma-glutamyl ligase":
        accession = "gamma-glutamyl ligase"
    if accession == "COG1538: Outer membrane protein":
        accession = "outer membrane protein"
    if accession == "FK-506 binding protein, peptidyl-prolyl cis-trans isomerase":
        accession = "peptidyl-prolyl cis-trans isomerase"
    if accession == "integral membrane protein SCI4115c":
        accession = "integral membrane protein"
    if accession == "G6PDH family F420-dependent oxidoreductase":
        accession = "oxidoreductase"

    accession = accession.replace("protein protein", "protein")

    # features contains 'golgi' but not contains 'golgi family' -> http://www.uniprot.org/uniprot/C5P901
    accession = accession.replace(
        "Golgi phosphoprotein 3", "Golgi phosphoprotein 3 family protein"
    )

    accession = accession.replace(", mitochondrial", "")
    if accession.find("mitochondrial") != -1:
        return ""

    # replace e. coli from the annotation
    accession = accession.replace("(E. coli)", "")

    if accession.lower().find("whole genome shotgun sequence") != -1:
        return ""

    accession = accession.strip()
    if accession == "Predicted protein":
        accession = "hypothetical protein"

    if accession.lower().startswith("regulator"):
        accession = "Regulator protein"

    """
        plural is not allowed
    """
    if accession.startswith("ATPase associated"):
        accession = "ATPase"
    elif accession == "transcriptional regulators":
        accession = "transcriptional regulator"

    if accession == "transferase":
        accession = "Transferase"
    accession = accession.strip(",")

    # mask tRNA/rRNA. so that i do not get removed
    accession = accession.replace("tRNA/rRNA", "tRNA / rRNA")

    # 'oxidoreductase/MSMEI_1564' -> 'oxidoreductase'
    if accession.find("dehydrogenase/reductase") > -1:
        accession = re.sub(
            "dehydrogenase/reductase/[a-zA-Z0-9_/.]+",
            "dehydrogenase/reductase",
            accession,
        )
    else:
        accession = re.sub("protein/[a-zA-Z0-9_/.]+", "protein", accession)
        accession = re.sub("regulator/[a-zA-Z0-9_/.]+", "regulator", accession)
        accession = re.sub("ter/[a-zA-Z0-9_/.]+", "ter", accession)
        accession = re.sub(
            "ase/[a-zA-Z0-9_/.]+", "ase", accession
        )  # like esterase/BLA654/FOO7888
        accession = re.sub("/[a-zA-Z0-9._]+/[a-zA-Z0-9/._]+", "", accession)

    # unmaks
    accession = accession.replace("tRNA / rRNA", "tRNA/rRNA")

    return accession


def change_glimmer3_prediction_output(infile, outfile):
    outfile = open(outfile, "w")
    contig_name = ""
    for line in open(infile):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            # Only take the Sequence ID
            contig_name = line[1:].split()[0]
            outfile.write(line + "\n")
        else:
            # orf-line, starts with the orfname at the beginning
            orf_name = line.split()[0]
            outfile.write(
                line.replace(orf_name, "%s_%s" % (contig_name, orf_name)) + "\n"
            )

    outfile.close()


if __name__ == "__main__":
    testdir = "./tests/"
    for filename in os.listdir(testdir):
        for line in open(os.path.join(testdir, filename)):
            if line.startswith("#"):
                continue
            in_, out_ = line.strip().split("\t")
            print(in_, out_, "-->", change_according_reviewer(in_.strip(), False))
            assert change_according_reviewer(in_.strip(), False) == out_
    sys.stdout.write("*** All tests passed. ***\n")
