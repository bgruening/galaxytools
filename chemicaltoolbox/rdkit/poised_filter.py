#!/usr/bin/env python

# Copyright 2017 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os

from rdkit import Chem
from rdkit.Chem import AllChem

from pipelines_utils_rdkit import rdkit_utils


class Filter(object):

    def __init__(self, *args, **kwargs):
        self.poised_filters = {
 'Amides': ['[#7:1][C;x0:2]=[O:3]>>[#7:1].Cl[C:2]=[O:3]'],
 'Benzimidazole': ['[#7:1]1[#6:9][#7:2][#6:4]2[#6:3]1[#6:8][#6:7][#6:6][#6:5]2>>[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#7:2].Cl[#6:9]=O'],
 'Benzoxazole': ['[#7:1]1[#6:9][#8:2][#6:4]2[#6:3]1[#6:8][#6:7][#6:6][#6:5]2>>[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#8:2].Cl[#6:9]=O'],
 'Ester_Coupling': ['[#8:1][C;x0:2]=[O:3]>>[#8:1].Cl[C:2]=[O:3]'],
 'Ether_Coupling': ['[CH0:2]-[#8R0:1]>>[#8:1].[#6:2]Br',
  '[CH1:2]-[#8R0:1]>>[#8:1].[#6:2]Br',
  '[CH2R0:2]-[#8:1]>>[#8:1].[#6:2]Br',
  '[CH1R0:2]-[#8:1]>>[#8:1].[#6:2]Br',
  '[CH2:2]-[#8R0:1]>>[#8:1].[#6:2]Br',
  '[CH3:2]-[#8R0:1]>>[#8:1].[#6:2]Br',
  '[CH0R0:2]-[#8:1]>>[#8:1].[#6:2]Br',
  '[n:3][c:2]-[#8:1]>>[#8:1].[n:3][#6:2]Br'],
 'Indole': ['[c:10]1[c:9][nH:1][c:3]2[c:8]1[c:7][c:6][c:5][c:4]2>>N[N:1][c:3]1[c:8]([2H])[c:7][c:6][c:5][c:4]1.[C:10][C:9]=O'],
 'N-Alkylation': ['[CH0:2]-[#7R0:1]>>[#7:1].[#6:2]Br',
  '[CH1:2]-[#7R0:1]>>[#7:1].[#6:2]Br',
  '[CH2R0:2]-[#7:1]>>[#7:1].[#6:2]Br',
  '[CH1R0:2]-[#7:1]>>[#7:1].[#6:2]Br',
  '[CH2:2]-[#7R0:1]>>[#7:1].[#6:2]Br',
  '[CH3:2]-[#7R0:1]>>[#7:1].[#6:2]Br',
  '[CH0R0:2]-[#7:1]>>[#7:1].[#6:2]Br'],
 'Oxadiazole': ['[#6:6][c:4]1[n:5][o:3][c:1][n:2]1>>[O:2]-[C:1]=[O:3].[#6:6][C:4]#[N:5]'],
 'Reductive_Amination': ['[CH1:2]-[#7R0:1]>>[#7:1].[#6:2]=O',
  '[CH2R0:2]-[#7:1]>>[#7:1].[#6:2]=O',
  '[CH1R0:2]-[#7:1]>>[#7:1].[#6:2]=O',
  '[CH2:2]-[#7R0:1]>>[#7:1].[#6:2]=O',
  '[CH3:2]-[#7R0:1]>>[#7:1].[#6:2]=O'],
 'SNAr': ['[c:2][N:1][#6:3]>>[#6:3]-[#7:1].[c:2]Br'],
 'Sonogashira': ['[#6;a:1][C:2]#[C:3]>>[#6;a:1]Br.[C:2]#[C:3]'],
 'Sulfonamide': ['[#7:1][S:2](=[O:3])=[O:4]>>[#7:1].Cl[S:2](=[O:3])=[O:4]'],
 'Suzuki_Coupling': ['[#6;a:1]-[#6;a:2]>>[#6;a:2]Br.[#6;a:1]-[#5](-[#8])-[#8]'],
 'Triazole': ['[#6:6][n:1]1[c:3][c:4][n:2][n:5]1>>[#6:6][#7:1][#7:5]=[#7:2].[C:3]#[C:4]'],
 'Urea': ['[#7:1][C;x0]([#7:2])=O>>[#7:1].[#7:2]']}

        # ++ means
        # THEN

        self.poised_reactions = {'Amides': ['[Cl,Br,I][C:2]=[O:3].[#7:1]>>[#7:1][C:2]=[O:3]',
                    '[OH1][C:2]=[O:3].[#7:1]>>[#7:1][C:2]=[O:3]',
                    '[NH2:1].Cl[C:2]=[O:3]>>[#7:1][C:2]=[O:3]',
                    '[NH1:1].Cl[C:2]=[O:3]>>[#7:1][C:2]=[O:3]'],
                                 'Sarah_Cu': [
                                     '[#6:1][N:2]=[N+:3]=[N-1:4].[C:5]#[C:6][C:7]>>[c:5]1=[c:6]([n:4]=[n:3][n:2]1[#6:1])[#6:7]',
                                           ],
                                 'Sarah_Quat_Am':[
                                     "[#6:1][C:2]([C:3](=[O:4])[O:5][C:6])([N:7])[#6:8].[c:9]1[c:10][c:11][c:12][c:13]([c:14]1[N:15])[C:16]([Cl,O:17])=[O:18]>>[c:9]1[c:10][c:11][c:12][c:13]2[c:14]1[N:15][C:3]([C:2]([N:7]([C:16]2=[O:18])[C:19])([#6:8])[#6:1])=[O:4]",
                                 ],
                                 #[c:1]1[c:2][c:3]([c:4]([c:5][c:6]1)[N:7])[C:8](=[O:9])[O,Cl:10]
                                 'SNAr': ['[#7H2:1].[c:2]Br>>[c:2][PH1:1]',
                  '[NH1:1].[c:2]Br>>[c:2][PH0:1]++[#15:1]>>[#7:1]',
                  '[c:2]-[Cl,Br,I].[#7:1]>>[c:2]-[#7:1]'],
                                 'Urea': ['[#7H2:1].[#7:2]>>[#7H1:1]C([#7:2])=O',
                  '[NH1:1].[#7:2]>>[#7H0:1]C([#7:2])=O',
                  '[#7:2]=C=O.[#7H2:1]>>[#7:1]-[#6](-[#7:2])=O',
                      '[#7:2]=C=O.[#7;AH1:1]>>[#7:1]-[#6](-[#7:2])=O'],
                                 'Suzuki Coupling':['[#6;a:2]-[#35,#53].[#6;a:1]-[#5](-[#8])-[#8]>>[#6;a:1]-[#6;a:2]',
        '[#6;a:1]-[#5](-[#8])-[#8].[#6;a:2]-[#35,#53]>>[#6;a:1]-[#6;a:2]'],
        'Sonogashira':['[C:1]#[C:2].[#6;a:3]-[#35,#53]>>[#6;a:3][C:1]#[C:2]',
        '[#6;a:1]-[#35,#53].[C:2]#[C:3]>>[#6;a:1][C:2]#[C:3]'],
                                 'Sulfonamide':['[#17,#35,#53][S:4](=[O:3])=[O:2].[#7:1]>>[#7:1][S:4](=[O:3])=[O:2]',
                                                '[#7H2:1].Cl[S:2](=[O:3])=[O:4]>>[#7:1][S:2](=[O:3])=[O:4]',
                                                '[#7H1:1].Cl[S:2](=[O:3])=[O:4]>>[#7:1][S:2](=[O:3])=[O:4]'],
                                 'Reductive_Amination':['[#6H1:2]=O.[#7:1]>>[#6:2]-[P:1]',
                                                        '[#6:3]-[#6:2](-[#6:4])=O.[#7:1]>>[#6:3]-[#6:2](-[#6:4])-[P:1]++[P:1]>>[N:1]',
                                                        '[#7;AH2:1].[#6:4]-[#6:2](-[#6:3])=O>>[#6:4]-[P:2](-[#6:3])-[#7:1]',
                                                        '[#7AH1:1].[#6:2]-[#6:3](-[#6:4])=O>>[#6:2]-[P:3](-[#6:4])-[#7:1]',
                                                        '[#7AH2:1].[#6:4]-[#6H1:2]=O>>[#6:4]-[P:2]-[#7:1]',
                                                        '[#7AH1:1].[#6:2]-[#6H1:3]=O>>[#6:2]-[P:3]-[#7:1]++[P:1]>>[C:1]'],
                                 'N-Alkylation':['[#7H1:1].[#6;A:2]Br>>[#7:1]-[#15:2]++[P:1]>>[C:1]',
                                                 '[C:2][Cl,Br,I].[#7H1:1]>>[#7:1]-[#15:2]++[P:1]>>[C:1]'],
                                 'Ether_Coupling':['[#8H1:1].[#6;A:2]Br>>[#6:2]-[#8:1]',
                                                   '[OH1:1].[n:3][#6;a:2]Br>>[O:1]-[#6:2][n:3]',
                                                   '[#6;A:2][#17,#35,#53].[#8H1:1]>>[#6:2]-[#8:1]',
                                                   '[n:3][c:2][#17,#35,#53].[#8H1:1]>>[n:3][c:2][#8:1]'],
                                 'Ester_Coupling':['[Cl,Br,I][C:2]=[O:3].[#8:1]>>[#8:1][C:2]=[O:3]',
                                                   '[OH1][C:2]=[O:3].[#8:1]>>[#8:1][C:2]=[O:3]',
                                                   '[OH:1].Cl[C:2]=[O:3]>>[#8:1][C:2]=[O:3]'],
                                 'Benzimidazole':['Cl[#6:9]=O.[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#7:2]>>[nH:1]1[c:9][n:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12',
                                                  '[OH][#6:9]=O.[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#7:2]>>[nH:1]1[c:9][n:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12',
                                                  '[N:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[N:2].Cl[#6:9]=O>>[nH:1]1[c:9][n:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12',
                                                  '[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#7:2].Cl[#6:9]=O>>[nH:1]1[c:9][n:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12'],
                                 'Triazole':['[C:3]#[C:4][#6:6].[#7:1][#7:5]=[#7:2]>>[#7:1]1[#6H:3]=[#6:4]([#6:6])[#7:2]=[#7:5]1',
                                             '[#7:1]=[#7+:2]=[#7-:3]>>[#7:1][#7:2]=[#7:3]++[#7:1][#7:5]=[#7:2].[C:3]#[C:4][#6:6]>>[#7:1]1[#6H:3]=[#6:4]([#6:6])[#7:2]=[#7:5]1'],
                                 'Benzoxazole':['Cl[#6:9]=O.[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#8:2]>>[n:1]1[c:9][o:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12',
                                                '[OH][#6:9]=O.[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#8:2]>>[n:1]1[c:9][o:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12',
                                                '[OH:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[NH2:2].Cl[#6:9]=O>>[o:1]1[c:9][n:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12']
        }
        self.reacts = {}
        self.starts = {}
        for key in self.poised_filters:
            self.starts[key] = [(Chem.MolFromSmarts(x.split(">>")[0]),AllChem.ReactionFromSmarts(x)) for x in self.poised_filters[key]]
        for key in self.poised_reactions:
            self.reacts[key] = []
            for x in self.poised_reactions[key]:
                react_seq = []
                for rxn in x.split("++"):
                    react_seq.append(AllChem.ReactionFromSmarts(rxn))
                self.reacts[key].append(react_seq)


    def pass_filter(self, mol):
        """
        Filter a given molecule given a series of SMARTS patterns.
        :param mol: an input RDKit Molecule
        :param patterns: a dict of SMARTS patterns
        :return: a list of the reactions each mol can do.
        """
        out_dict = {}
        for key in self.starts:
            for patt_rxn_pair in self.starts[key]:
                patt = patt_rxn_pair[0]
                rxn = patt_rxn_pair[1]
                if mol.HasSubstructMatch(patt):
                    products = self.unique_products(rxn.RunReactants([mol]))
                    if key in out_dict:
                        out_dict[key].extend(products)
                        out_dict[key] = list(set(out_dict[key]))
                    else:
                        out_dict[key] = products
        return out_dict


    def unique_products(self, ps, min_num_heavy_atoms=2):
        """
        Get all the unique products out
        :param ps: the input reaction products
        :param min_num_heavy_atoms: the minimum number of heavy atoms
        :return:
        """
        uniqps = {}
        for products in ps:
            for p in products:
                #### Filter out singletons
                if p.GetNumHeavyAtoms() >= min_num_heavy_atoms:
                    smi = Chem.MolToSmiles(p, isomericSmiles=True)
                    uniqps[smi] = p
        products = sorted(uniqps.keys())
        return products



    def get_writers(self, dir_base):
        """
        Get all the writers of the SD files
        :param output_path:
        :return:
        """
        out_d = {}
        for x in self.starts:
            out_d[x] = Chem.SDWriter(os.path.join(dir_base, x + ".sdf"))
        return out_d

    def annotate_reagent(self, product, reactants, reaction_name):
        """
        Annotate the product with the input molecules (and their input molecules).
        :param product:
        :param reactants:
        :param reaction_name:
        :return:
        """
        reagent_delimiter = "WITH"
        reaction_namer = "BYE"
        reaction_delimiter = "AND"

        smis = [Chem.MolToSmiles(x,isomericSmiles=True) for x in reactants]
        current_reaction = reagent_delimiter.join(smis)+reaction_namer+reaction_name
        provenances = [x.GetProp("SOURCE") for x in reactants if x.HasProp("SOURCE")]
        product.SetProp("SOURCE",reaction_delimiter.join(provenances)+reaction_delimiter+current_reaction)



    def run_reaction(self,input_molecule,reactant_mol, react_seqs):
        products = []
        for react_seq in react_seqs:
            rxn = react_seq[0]
            products.extend(rxn.RunReactants((input_molecule, reactant_mol,)))
            for i in range(1,len(react_seq)):
                rxn = react_seq[i]
                curr_prods = products
                products = []
                for p in curr_prods:
                    products.extend(rxn.RunReactants(p,))
        return products

    def perform_reaction(self, input_molecule, reaction_name, reactant_mol, writer,i):
        """Take an input molecule and a library of reactants
        React - and form the products.
        :param input_molecule: the input molecule to be reacted
        :param reaction_name: the name of the reaction (a string)
        :param reactant_lib: an iterable of molecules to react
        :return:
        """
        if input_molecule.HasProp('uuid'):
            mol_uuid = input_molecule.GetProp('uuid')
        else:
            mol_uuid = None
        react_seq = self.reacts[reaction_name]
        if reactant_mol.HasProp('uuid'):
            reactant_uuid = reactant_mol.GetProp('uuid')
        else:
            reactant_uuid = None
        products = [Chem.MolFromSmiles(x) for x in self.unique_products(self.run_reaction(input_molecule,reactant_mol,react_seq)) if Chem.MolFromSmiles(x)]
        for product in products:
            i+=1
            rdkit_utils.generate_2d_coords(product)
            if mol_uuid:
                product.SetProp("source_uuid", mol_uuid)
            if reactant_uuid:
                product.SetProp("reactant_uuid", reactant_uuid)
            writer.write(product)
        return i

    def get_subs(self,mol,reaction_name):
        """
        Get all the substructures for a given reaction.
        :param mol:
        :param reaction_name:
        :return:
        """
        return mol.GetProp(reaction_name).split(",")

    def get_rxn_names(self):
        return self.poised_filters.keys()