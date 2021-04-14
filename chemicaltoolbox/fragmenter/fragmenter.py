#!/usr/bin/env python
import openbabel

openbabel.obErrorLog.StopLogging()
import pybel
import sys
import argparse

"""
Based on a script from TJ O'Donnell
https://gist.github.com/95387

Modified by B. Gruening and Hitesh Patel
Copyright 2009-2012 TJ O'Donnell, B.Gruening
"""


class Fragmenter:
    def __init__(self, mol, minsize=5):
        self.mol = mol
        # minimum allowed size (atom count) of fragment
        self.minsize = minsize
        # bonded atom pairs populated by the apply method,
        # subsequently used by split and add_star
        self.atom_pairs = list()
        self.parent_name = mol.title or ""
        self.parent_can = mol.write("can").strip() or ""

    def __set_atomic_num(
        self, atom1, atom2, first_mark, second_mark, check_first_atomic_num=False
    ):
        """
        Change the atomic number of the sticky ends to mark it for further calculations.
        first_mark = The atomic number that should mark the first sticky end that corresponds to the check_atomic_num.
        second_mark = The atomic number that should mark the second sticky end.
        check_atomic_num = The original atomic number for the first sticky end.
        """
        if check_first_atomic_num and check_first_atomic_num == atom1.GetAtomicNum():
            atom1.SetAtomicNum(first_mark)
            atom2.SetAtomicNum(second_mark)
        else:
            atom1.SetAtomicNum(second_mark)
            atom2.SetAtomicNum(first_mark)

    def __mark_atoms(self, pattern_num, atom1, atom2):
        """
        Replaces the atoms the new end-atoms that are left from the bond
        break with  some really heavy atoms (atomic number from 89-108)
        to mark them. In furhter steps these marks are recognised and a new
        bond is formed.
        """
        if pattern_num == 1:
            # 1C 1N
            self.__set_atomic_num(atom1, atom2, 89, 90, 6)
        elif pattern_num == 2:
            # 2C 2O
            self.__set_atomic_num(atom1, atom2, 91, 92, 6)
        elif pattern_num == 3:
            # 3C 3N
            self.__set_atomic_num(atom1, atom2, 93, 94, 6)
        elif pattern_num == 4:
            # 4C 4N
            self.__set_atomic_num(atom1, atom2, 95, 96, 6)
        elif pattern_num == 5:
            # 5C 5O
            self.__set_atomic_num(atom1, atom2, 97, 98, 6)
        elif pattern_num == 6:
            # 6C 6C
            self.__set_atomic_num(atom1, atom2, 99, 99)
        elif pattern_num == 7:
            # 7C 7N
            self.__set_atomic_num(atom1, atom2, 100, 101, 6)
        elif pattern_num == 8:
            # 8C 8N
            self.__set_atomic_num(atom1, atom2, 102, 103, 6)
        elif pattern_num == 9:
            # 9C 9N
            self.__set_atomic_num(atom1, atom2, 104, 105, 6)
        elif pattern_num == 10:
            # 10C 10C
            self.__set_atomic_num(atom1, atom2, 106, 106)
        elif pattern_num == 11:
            # 11N 11S
            self.__set_atomic_num(atom1, atom2, 107, 108, 7)

    def fragmentation(self, pattern, pattern_num):
        hits = pattern.findall(self.mol)
        if hits:
            # find all atom pairs that match
            for atom_id_one, atom_id_two in hits:
                if atom_id_one < atom_id_two:
                    atom_one = self.mol.atoms[atom_id_one - 1].OBAtom
                    atom_two = self.mol.atoms[atom_id_two - 1].OBAtom
                else:
                    atom_two = self.mol.atoms[atom_id_one - 1].OBAtom
                    atom_one = self.mol.atoms[atom_id_two - 1].OBAtom

                if not self.is_fragment_to_small(atom_one, atom_two):
                    self.atom_pairs.append([atom_one, atom_two, pattern_num])

    def split(self, mark_sticky_ends):
        for atom_one, atom_two, pattern_num in self.atom_pairs:
            bond = atom_one.GetBond(atom_two)
            if bond:
                if mark_sticky_ends:
                    self.__mark_atoms(pattern_num, atom_one, atom_two)

                self.mol.OBMol.DeleteBond(bond)
            else:
                sys.stderr.write(
                    "Could not split %s at bond (%s, %s)\n"
                    % (
                        self.mol.write("can").split()[0],
                        atom_one.GetIdx(),
                        atom_two.GetIdx(),
                    )
                )
                return False

    def decide_multiples(self):
        """
        some smarts (e.g. ether, amine) allow multiple bonds to the
        central atom to be broken. Yet it appears the central atom
        needs to be retained in one of the multiple fragments.
        If multiple fragments, let it stay with the smallest fragment.
        If tied, pick the first fragment.
        """
        # TODO: check why the +1 is necessary
        multiples = [0] * (self.mol.OBMol.NumAtoms() + 1)
        for pair in self.atom_pairs:
            multiples[pair[0].GetIdx()] += 1
            multiples[pair[1].GetIdx()] += 1

        while max(multiples) > 1:
            """
            We need to itereate over that construct more than one time, when we match one atom more than two times.
            For example: S(=O)(=O)(Nc1onc(c1C)C)c1ccc(N)cc1
            multiples are -> [0, 2, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            """
            currsize = -1
            currpair = None
            for pair in self.atom_pairs:
                a, b, pattern_num = pair
                if multiples[a.GetIdx()] > 1 or multiples[b.GetIdx()] > 1:
                    # remove larger fragment(s) if a-b were broken
                    fsize = self.fragment_size(a, b)
                    if currpair == None:
                        currpair = pair
                        currsize = fsize
                    else:
                        if fsize < currsize:
                            self.atom_pairs.remove(pair)
                            multiples[a.GetIdx()] -= 1
                            multiples[b.GetIdx()] -= 1
                        else:
                            self.atom_pairs.remove(currpair)
                            multiples[currpair[0].GetIdx()] -= 1
                            multiples[currpair[1].GetIdx()] -= 1
                            currpair = pair
                            currsize = fsize

    def fragment_size(self, a, b):
        # size of fragment b if a-b were broken
        c1 = openbabel.vectorInt()
        self.mol.OBMol.FindChildren(c1, a.GetIdx(), b.GetIdx())
        return len(c1) + 1

    def is_fragment_to_small(self, a, b):
        # if we were to break the bond between a and b,
        #    would either fragment be too small?
        if self.fragment_size(a, b) < self.minsize:
            return True
        if self.fragment_size(b, a) < self.minsize:
            return True

        return False

    def separate(self, molecule_type, outfile):
        output = pybel.Outputfile(molecule_type, outfile, overwrite=True)
        split_mol = self.mol.OBMol.Separate()

        for mol in split_mol:
            mol = pybel.Molecule(mol)
            mol.data.update(
                {
                    "parent name": self.parent_name,
                    "parent canonical smiles": self.parent_can,
                }
            )
            output.write(mol)
        output.close()

    def get_fragments(self):
        fragments = []
        for mol in self.mol.OBMol.Separate():
            mol = pybel.Molecule(mol)
            mol.data.update(
                {
                    "parent name": self.parent_name,
                    "parent canonical smiles": self.parent_can,
                }
            )
            fragments.append(mol)
        return fragments


def __main__():
    parser = argparse.ArgumentParser(
        description="Splits a molecule to several fragments."
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="input_path",
        required=True,
        help="Path to the input file.",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output_path",
        required=True,
        help="Path to the output file.",
    )

    parser.add_argument(
        "-n",
        "--non-fragments",
        dest="non_fragment_file",
        help="Specify an additional file to store all molecules that aren't fragmented at all. Without that option such molecules will also appear in the regular output file.",
    )

    parser.add_argument(
        "--input-format",
        dest="iformat",
        help="Input format. It must be supported by openbabel.",
    )

    parser.add_argument(
        "--output-format",
        dest="oformat",
        default="can",
        help="Output format. It must be supported by openbabel.",
    )

    parser.add_argument(
        "-m",
        "--mark-sticky-ends",
        dest="mark_sticky_ends",
        action="store_true",
        default=False,
        help="Replaces the newly created free ends on each bond break with some marker atoms. "
        "This option is usefull if you want to to merge the fragments afterwarts.",
    )

    parser.add_argument(
        "-r",
        "--splitting-rules",
        dest="rules",
        choices=["recap", "rotbonds", "rings"],
        default="recap",
        help="Splitting rules.",
    )

    # each smarts must contain only two atoms representing the bond to be broken.
    # of course, each atom may be a complex atom smarts, ala [$(whatever)]
    """
        RECAP-Retrosynthetic Combinatorial Analysis Procedure
        J. Chem. Inf. Comput. Sci. 1998, 38, 511-522
    """
    recap_smarts = {
        " 1.amide": "[$([C;!$(C([#7])[#7])](=!@[O]))]!@[$([#7;+0;!D1])]",
        " 2.ester": "[$(C=!@O)]-!@[$([O;+0])]",
        " 3.amine": "[$([N;!D1;+0;!$(N-C=[#7,#8,#15,#16])](-!@[*]))]-!@[$([#6])]",
        " 4.urea": "[$(C(=!@O)([#7;+0;D2,D3])!@[#7;+0;D2,D3])]!@[$([#7;+0;D2,D3])]",
        " 5.ether": "[$([O;+0](-!@[#6!$(C=O)])-!@[#6!$(C=O)])]-!@[$([#6!$(C=O)])]",
        " 6.olefin": "C=!@C",
        " 7.quaternaryN": "[N;+1;D4]!@[#6]",
        " 8.aromaticN-aliphaticC": "[$([n;+0])]-!@C",
        " 9.lactamN-aromaticC": "[$([N;+0]-@[C]=[O])]-!@[$([c])]",
        "10.aromaticC-aromaticC": "c-!@c",
        "11.sulphonamide": "[$([#7;+0;D2,D3])]-!@[$([S](=[O])=[O])]",
    }

    """
    recap_smarts = {
    ' 1.amide':'[$([C;!$(C([#7])[#7])](=!@[O]))]!@[$([#7;+0;!D1])]',
    ' 2.ester':'[$(C=!@O)]!@[$([O;+0])]',
    ' 3.amine':'[$([N;!D1;+0;!$(N-C=[#7,#8,#15,#16])](-!@[*]))]-!@[$([*])]',
    ' 4.urea':'[$(C(=!@O)([#7;+0;D2,D3])!@[#7;+0;D2,D3])]!@[$([#7;+0;D2,D3])]',
    ' 5.ether':'[$([O;+0](-!@[#6!$(C=O)])-!@[#6!$(C=O)])]-!@[$([#6!$(C=O)])]',
    ' 6.olefin':'C=!@C',
    ' 7.quaternaryN':'[N;+1;D4]!@[#6]',
    ' 8.aromaticN-aliphaticC':'[$([n;+0])]-!@C',
    ' 9.lactamN-aromaticC':'[$([O]=[C]-@[N;+0])]-!@[$([C])]',
    '10.aromaticC-aromaticC':'c-!@c',
    '11.sulphonamide':'[$([#7;+0;D2,D3])]-!@[$([S](=[O])=[O])]',
    }

    """
    rotatable_bonds_smarts = {"1.rotatable-bonds": "[!$(*#*)&!D1]-!@[!$(*#*)&!D1]"}

    ring_smarts = {
        "1.non-ring-bond-between-to-rings": "[R]!@[R]",
        "2.ring-associated-atom": "[R]~[!R]",
    }

    options = parser.parse_args()
    output = pybel.Outputfile(options.oformat, options.output_path, overwrite=True)

    if options.non_fragment_file:
        non_fragment_file = pybel.Outputfile(
            options.oformat, options.non_fragment_file, overwrite=True
        )

    if not options.iformat:
        from cheminfolib import check_filetype

        options.iformat = check_filetype(options.input_path)

    for mol in pybel.readfile(options.iformat, options.input_path):
        frag = Fragmenter(mol, 4)
        if options.rules == "rotbonds":
            smarts = rotatable_bonds_smarts
        elif options.rules == "rings":
            smarts = ring_smarts
        else:
            smarts = recap_smarts
        for i, stype in enumerate(sorted(smarts.keys()), start=1):
            frag.fragmentation(pybel.Smarts(smarts[stype]), i)
        frag.decide_multiples()
        frag.split(options.mark_sticky_ends)

        fragments = frag.get_fragments()

        # if a separate file is specified for molecules without a bond-break,
        # write it to that file and continue with the next molecule
        if options.non_fragment_file and len(fragments) == 1:
            non_fragment_file.write(fragments[0])
            continue

        for mol in fragments:
            # write in a given molecule format
            output.write(mol)
    output.close()


if __name__ == "__main__":
    __main__()
