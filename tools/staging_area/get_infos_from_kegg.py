#!/usr/bin/env python
"""
Input: DNA Fasta File
Output: Tabular
Return Tabular File with predicted ORF's
Bjoern Gruening
"""
import sys, os
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.expression import func, cast, alias
from kegg import KEGGDB


KEGGDBEngine, KEGGDBBase = KEGGDB.init('postgres://bag:bag@10.4.56.8:5432/phabidb')
KEGGSession = sessionmaker(bind=KEGGDBEngine)()


def __main__():

    infile_path = sys.argv[1]
    column = sys.argv[2]
    outfile_path = sys.argv[3]

    outfile = open(outfile_path, 'w')

    for line in open(infile_path, 'r'):
        token = line.split('\t')
        gene_id = token[int(column) -1].strip()
        results = KEGGSession.query(KEGGDB.GenesProteinSeq).filter(
                KEGGDB.GenesProteinSeq.gene_id.op('LIKE')('%'+gene_id+'%')
            ).all()

        data = []
        for result in results:
            data.append( ' '.join([ko.ko for ko in result.kos]) )
            data.append( ' '.join([ec.ec for ec in result.ecs]) )
            data.append( ' '.join([pathway.pathway for pathway in result.pathways]) )
            data.append( ' '.join([uniprot.uniprot for uniprot in result.uniprots]) )
            data.append( ' '.join([prosite.prosite for prosite in result.prosites]) )
            data.append( ' '.join([ncbi.ncbi_geneid for ncbi in result.ncbigenes]) )
            data.append( ' '.join([pfam.pfam for pfam in result.pfams]) )
        outfile.write(line.strip() + "\t" + '\t'.join(data)+'\n')

if __name__ == "__main__" :
    __main__()
