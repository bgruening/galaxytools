wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/gp2protein/gp2protein.unigene.gz
gunzip gp2protein.unigene.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/text/pdb_chain_pfam.lst

wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
gunzip idmapping.dat.gz

wget http://reactome.oicr.on.ca/download/current/uniprot_2_pathways.stid.txt

wget ftp://ftp.ncbi.nih.gov/pubchem/Compound/Extras/CID-Parent.gz
gunzip CID-Parent.gz

wget ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Synonym-unfiltered.gz
gunzip CID-Synonym-unfiltered.gz

wget ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Synonym-filtered.gz
gunzip CID-Synonym-filtered.gz

wget ftp://ftp.ncbi.nih.gov/pubchem/Compound/Extras/CID-MeSH


wget ftp://ftp.ncbi.nih.gov/pubchem/Compound/Extras/MeSH-Pharm

#This file contains all GO annotations for proteins in the UniProt KnowledgeBase (UniProtKB).
wget ftp://ftp.geneontology.org/pub/go/gene-associations/submission/gene_association.goa_uniprot.gz
gunzip gene_association.goa_uniprot.gz

wget http://ligand-expo.rcsb.org/dictionaries/cc-to-pdb.tdd

wget http://research.isb-sib.ch/unimed/SP_MeSH.tab

wget http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.1.85/BIOGRID-IDENTIFIERS-3.1.85.tab.zip
unzip BIOGRID-IDENTIFIERS-3.1.85.tab.zip

wget http://research.isb-sib.ch/unimed/SP_MeSH.tab


cd go

wget http://www.geneontology.org/external2go/ec2go
wget http://www.geneontology.org/external2go/interpro2go
wget http://www.geneontology.org/external2go/kegg2go
wget http://www.geneontology.org/external2go/pfam2go
wget http://www.geneontology.org/external2go/reactome2go
wget http://www.geneontology.org/doc/GO.terms_alt_ids

cd ..
