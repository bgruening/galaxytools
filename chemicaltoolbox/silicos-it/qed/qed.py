#!/usr/bin/env python
__all__ = ['weights_max', 'weights_mean', 'weights_none', 'default']

# Silicos-it
from errors import WrongArgument

# RDKit
from rdkit.Chem import Descriptors
from rdkit import Chem

# General
from copy import deepcopy
from math import exp, log
import sys, os, re
import argparse

def check_filetype(filepath):
    mol = False
    possible_inchi = True
    for line_counter, line in enumerate(open(filepath)):
        if line_counter > 10000:
            break
        if line.find('$$$$') != -1:
            return 'sdf'
        elif line.find('@<TRIPOS>MOLECULE') != -1:
            return 'mol2'
        elif line.find('ligand id') != -1:
            return 'drf'
        elif possible_inchi and re.findall('^InChI=', line):
            return 'inchi'
        elif re.findall('^M\s+END', line):
            mol = True
        # first line is not an InChI, so it can't be an InChI file
        possible_inchi = False

    if mol:
        # END can occures before $$$$, so and SDF file will 
        # be recognised as mol, if you not using this hack'
        return 'mol'
    return 'smi'

AliphaticRings = Chem.MolFromSmarts('[$([A;R][!a])]')

AcceptorSmarts = [
    '[oH0;X2]',
    '[OH1;X2;v2]',
    '[OH0;X2;v2]',
    '[OH0;X1;v2]',
    '[O-;X1]',
    '[SH0;X2;v2]',
    '[SH0;X1;v2]',
    '[S-;X1]',
    '[nH0;X2]',
    '[NH0;X1;v3]',
    '[$([N;+0;X3;v3]);!$(N[C,S]=O)]'
    ]
Acceptors = []
for hba in AcceptorSmarts:
    Acceptors.append(Chem.MolFromSmarts(hba))

StructuralAlertSmarts = [
    '*1[O,S,N]*1',
    '[S,C](=[O,S])[F,Br,Cl,I]',
    '[CX4][Cl,Br,I]',
    '[C,c]S(=O)(=O)O[C,c]',
    '[$([CH]),$(CC)]#CC(=O)[C,c]',
    '[$([CH]),$(CC)]#CC(=O)O[C,c]',
    'n[OH]',
    '[$([CH]),$(CC)]#CS(=O)(=O)[C,c]',
    'C=C(C=O)C=O',
    'n1c([F,Cl,Br,I])cccc1',
    '[CH1](=O)',
    '[O,o][O,o]',
    '[C;!R]=[N;!R]',
    '[N!R]=[N!R]',
    '[#6](=O)[#6](=O)',
    '[S,s][S,s]',
    '[N,n][NH2]',
    'C(=O)N[NH2]',
    '[C,c]=S',
    '[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]',
    'C1(=[O,N])C=CC(=[O,N])C=C1',
    'C1(=[O,N])C(=[O,N])C=CC=C1',
    'a21aa3a(aa1aaaa2)aaaa3',
    'a31a(a2a(aa1)aaaa2)aaaa3',
    'a1aa2a3a(a1)A=AA=A3=AA=A2',
    'c1cc([NH2])ccc1',
    '[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]',
    'I',
    'OS(=O)(=O)[O-]',
    '[N+](=O)[O-]',
    'C(=O)N[OH]',
    'C1NC(=O)NC(=O)1',
    '[SH]',
    '[S-]',
    'c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]',
    'c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]',
    '[CR1]1[CR1][CR1][CR1][CR1][CR1][CR1]1',
    '[CR1]1[CR1][CR1]cc[CR1][CR1]1',
    '[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1',
    '[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1',
    '[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1',
    '[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1',
    'C#C',
    '[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]',
    '[$([N+R]),$([n+R]),$([N+]=C)][O-]',
    '[C,c]=N[OH]',
    '[C,c]=NOC=O',
    '[C,c](=O)[CX4,CR0X3,O][C,c](=O)',
    'c1ccc2c(c1)ccc(=O)o2',
    '[O+,o+,S+,s+]',
    'N=C=O',
    '[NX3,NX4][F,Cl,Br,I]',
    'c1ccccc1OC(=O)[#6]',
    '[CR0]=[CR0][CR0]=[CR0]',
    '[C+,c+,C-,c-]',
    'N=[N+]=[N-]',
    'C12C(NC(N1)=O)CSC2',
    'c1c([OH])c([OH,NH2,NH])ccc1',
    'P',
    '[N,O,S]C#N',
    'C=C=O',
    '[Si][F,Cl,Br,I]',
    '[SX2]O',
    '[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)',
    'O1CCCCC1OC2CCC3CCCCC3C2',
    'N=[CR0][N,n,O,S]',
    '[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2]2',
    'C=[C!r]C#N',
    '[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1',
    '[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1',
    '[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])',
    '[OH]c1ccc([OH,NH2,NH])cc1',
    'c1ccccc1OC(=O)O',
    '[SX2H0][N]',
    'c12ccccc1(SC(S)=N2)',
    'c12ccccc1(SC(=S)N2)',
    'c1nnnn1C=O',
    's1c(S)nnc1NC=O',
    'S1C=CSC1=S',
    'C(=O)Onnn',
    'OS(=O)(=O)C(F)(F)F',
    'N#CC[OH]',
    'N#CC(=O)',
    'S(=O)(=O)C#N',
    'N[CH2]C#N',
    'C1(=O)NCC1',
    'S(=O)(=O)[O-,OH]',
    'NC[F,Cl,Br,I]',
    'C=[C!r]O',
    '[NX2+0]=[O+0]',
    '[OR0,NR0][OR0,NR0]',
    'C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]',
    '[CX2R0][NX3R0]',
    'c1ccccc1[C;!R]=[C;!R]c2ccccc2',
    '[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]',
    '[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]',
    '[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]',
    '[*]=[N+]=[*]',
    '[SX3](=O)[O-,OH]',
    'N#N',
    'F.F.F.F',
    '[R0;D2][R0;D2][R0;D2][R0;D2]',
    '[cR,CR]~C(=O)NC(=O)~[cR,CR]',
    'C=!@CC=[O,S]',
    '[#6,#8,#16][C,c](=O)O[C,c]',
    'c[C;R0](=[O,S])[C,c]',
    'c[SX2][C;!R]',
    'C=C=C',
    'c1nc([F,Cl,Br,I,S])ncc1',
    'c1ncnc([F,Cl,Br,I,S])c1',
    'c1nc(c2c(n1)nc(n2)[F,Cl,Br,I])',
    '[C,c]S(=O)(=O)c1ccc(cc1)F',
    '[15N]',
    '[13C]',
    '[18O]',
    '[34S]'
    ]

StructuralAlerts = []
for smarts in StructuralAlertSmarts:
    StructuralAlerts.append(Chem.MolFromSmarts(smarts))


# ADS parameters for the 8 molecular properties: [row][column]
#     rows[8]:     MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS
#     columns[7]: A, B, C, D, E, F, DMAX
# ALOGP parameters from Gregory Gerebtzoff (2012, Roche)
pads1 = [    [2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677, 65.37051707, 104.9805561],
            [0.486849448, 186.2293718, 2.066177165, 3.902720615, 1.027025453, 0.913012565, 145.4314800],
            [2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953, 1.300669958, 148.7763046],    
            [1.618662227, 1010.051101, 0.985094388, 0.000000001, 0.713820843, 0.920922555, 258.1632616],
            [1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824, 28.51324732, 104.5686167],
            [0.010000000, 272.4121427, 2.558379970, 1.565547684, 1.271567166, 2.758063707, 105.4420403],
            [3.217788970, 957.7374108, 2.274627939, 0.000000001, 1.317690384, 0.375760881, 312.3372610],
            [0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477, 0.875193782, 417.7253140]        ]
# ALOGP parameters from the original publication
pads2 = [    [2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677, 65.37051707, 104.9805561],
            [3.172690585, 137.8624751, 2.534937431, 4.581497897, 0.822739154, 0.576295591, 131.3186604],
            [2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953, 1.300669958, 148.7763046],    
            [1.618662227, 1010.051101, 0.985094388, 0.000000001, 0.713820843, 0.920922555, 258.1632616],
            [1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824, 28.51324732, 104.5686167],
            [0.010000000, 272.4121427, 2.558379970, 1.565547684, 1.271567166, 2.758063707, 105.4420403],
            [3.217788970, 957.7374108, 2.274627939, 0.000000001, 1.317690384, 0.375760881, 312.3372610],
            [0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477, 0.875193782, 417.7253140]        ]

def ads(x, a, b, c, d, e, f, dmax):
    return ((a+(b/(1+exp(-1*(x-c+d/2)/e))*(1-1/(1+exp(-1*(x-c-d/2)/f))))) / dmax)

def properties(mol):
    """
    Calculates the properties that are required to calculate the QED descriptor.
    """
    matches = []
    if mol is None:
        raise WrongArgument("properties(mol)", "mol argument is \'None\'")
    x = [0] * 9
    x[0] = Descriptors.MolWt(mol)                                                # MW 
    x[1] = Descriptors.MolLogP(mol)                                                # ALOGP
    for hba in Acceptors:                                                        # HBA
        if mol.HasSubstructMatch(hba):
            matches = mol.GetSubstructMatches(hba)
            x[2] += len(matches)
    x[3] = Descriptors.NumHDonors(mol)                                            # HBD
    x[4] = Descriptors.TPSA(mol)                                                # PSA
    x[5] = Descriptors.NumRotatableBonds(mol)                                    # ROTB
    x[6] = Chem.GetSSSR(Chem.DeleteSubstructs(deepcopy(mol), AliphaticRings))    # AROM
    for alert in StructuralAlerts:                                                # ALERTS
        if (mol.HasSubstructMatch(alert)): x[7] += 1
    ro5_failed = 0
    if x[3] > 5:
        ro5_failed += 1 #HBD
    if x[2] > 10:
        ro5_failed += 1 #HBA
    if x[0] >= 500:
        ro5_failed += 1
    if x[1] > 5:
        ro5_failed += 1
    x[8] = ro5_failed
    return x


def qed(w, p, gerebtzoff):
    d = [0.00] * 8
    if gerebtzoff:
        for i in range(0, 8):
            d[i] = ads(p[i], pads1[i][0], pads1[i][1], pads1[i][2], pads1[i][3], pads1[i][4], pads1[i][5], pads1[i][6])
    else:
        for i in range(0, 8):
            d[i] = ads(p[i], pads2[i][0], pads2[i][1], pads2[i][2], pads2[i][3], pads2[i][4], pads2[i][5], pads2[i][6])
    t = 0.0
    for i in range(0, 8):
        t += w[i] * log(d[i])
    return (exp(t / sum(w)))


def weights_max(mol, gerebtzoff = True, props = False):
    """
    Calculates the QED descriptor using maximal descriptor weights.
    If props is specified we skip the calculation step and use the props-list of properties.
    """
    if not props:
        props = properties(mol)
    return qed([0.50, 0.25, 0.00, 0.50, 0.00, 0.50, 0.25, 1.00], props, gerebtzoff)


def weights_mean(mol, gerebtzoff = True, props = False):
    """
    Calculates the QED descriptor using average descriptor weights.
    If props is specified we skip the calculation step and use the props-list of properties.
    """
    if not props:
        props = properties(mol)
    return qed([0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95], props, gerebtzoff)


def weights_none(mol, gerebtzoff = True, props = False):
    """
    Calculates the QED descriptor using unit weights.
    If props is specified we skip the calculation step and use the props-list of properties.
    """
    if not props:
        props = properties(mol)
    return qed([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00], props, gerebtzoff)


def default(mol, gerebtzoff = True):
    """
    Calculates the QED descriptor using average descriptor weights and Gregory Gerebtzoff parameters.
    """
    return weights_mean(mol, gerebtzoff)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', 
                required=True, 
                help='path to the input file name')
    parser.add_argument("-m", "--method", 
                dest="method",
                choices=['max', 'mean', 'unweighted'],
                default="mean",
                help="Specify the method you want to use.")
    parser.add_argument("--iformat",
                help="Input format. It must be supported by openbabel.")
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w+'), 
                default=sys.stdout, 
                help="path to the result file, default it sdtout")
    parser.add_argument("--header", dest="header", action="store_true",
                default=False,
                help="Write header line.")


    args = parser.parse_args()

    # Elucidate filetype and open supplier
    ifile = os.path.abspath(args.input)
    if not os.path.isfile(ifile):
        print "Error: ", ifile, " is not a file or cannot be found."
        sys.exit(1)
    if not os.path.exists(ifile):
        print "Error: ", ifile, " does not exist or cannot be found."
        sys.exit(1)
    if not os.access(ifile, os.R_OK):
        print "Error: ", ifile, " is not readable."
        sys.exit(1)

    if not args.iformat:
        # try to guess the filetype
        filetype = check_filetype( ifile )
    else:
        filetype = args.iformat # sdf or smi


    """
        We want to store the original SMILES in the output. So in case of a SMILES file iterate over the file and convert each line separate.
    """
    if filetype == 'sdf':
        supplier = Chem.SDMolSupplier( ifile )
        # Process file
        if args.header:
            args.outfile.write("MW\tALOGP\tHBA\tHBD\tPSA\tROTB\tAROM\tALERTS\tLRo5\tQED\tNAME\n")
        count = 0
        for mol in supplier:
            count += 1
            if mol is None:
                print "Warning: skipping molecule ", count, " and continuing with next."
                continue
            props = properties(mol)

            if args.method == 'max':
                calc_qed = weights_max(mol, True, props)
            elif args.method == 'unweighted':
                calc_qed = weights_none(mol, True, props)
            else:
                calc_qed = weights_mean(mol, True, props)

            args.outfile.write( "%.2f\t%.3f\t%d\t%d\t%.2f\t%d\t%d\t%d\t%s\t%.3f\t%-s\n" % (
                props[0], 
                props[1], 
                props[2], 
                props[3], 
                props[4], 
                props[5], 
                props[6], 
                props[7],
                props[8],
                calc_qed,
                mol.GetProp("_Name"),
                ))
    elif filetype == 'smi':
        supplier = Chem.SmilesMolSupplier( ifile, " \t", 0, 1, False, True )

        # Process file
        if args.header:
            args.outfile.write("MW\tALOGP\tHBA\tHBD\tPSA\tROTB\tAROM\tALERTS\tLRo5\tQED\tNAME\tSMILES\n")
        count = 0
        for line in open(ifile):
            tokens = line.strip().split('\t')
            if len(tokens) > 1:
                smiles, title = tokens
            else:
                smiles = tokens[0]
                title = ''
            mol = Chem.MolFromSmiles(smiles)
            count += 1
            if mol is None:
                print "Warning: skipping molecule ", count, " and continuing with next."
                continue
            props = properties(mol)

            if args.method == 'max':
                calc_qed = weights_max(mol, True, props)
            elif args.method == 'unweighted':
                calc_qed = weights_none(mol, True, props)
            else:
                calc_qed = weights_mean(mol, True, props)

            args.outfile.write( "%.2f\t%.3f\t%d\t%d\t%.2f\t%d\t%d\t%d\t%s\t%.3f\t%-s\t%s\n" % (
                props[0], 
                props[1], 
                props[2], 
                props[3], 
                props[4], 
                props[5], 
                props[6], 
                props[7],
                props[8],
                calc_qed,
                title,
                smiles
                ))
    else:
        sys.exit("Error: unknown file-type: %s" % filetype)
