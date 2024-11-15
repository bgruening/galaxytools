## OSRA: Optical Structure Recognition Application

OSRA is a utility designed to convert graphical representations of chemical 
structures, as they appear in journal articles, patent documents, textbooks, 
trade magazines etc., into SMILES (Simplified Molecular Input Line Entry 
Specification - see http://en.wikipedia.org/wiki/SMILES) or 
SD files - a computer recognizable molecular structure format. 
OSRA can read a document in any of the over 90 graphical formats parseable by 
ImageMagick - including GIF, JPEG, PNG, TIFF, PDF, PS etc., and generate 
the SMILES or SDF representation of the molecular structure images encountered 
within that document.

Note that any software designed for optical recognition is unlikely to be 
perfect, and the output produced might, and probably will, contain errors, 
so curation by a human knowledgeable in chemical structures is highly recommended.

http://cactus.nci.nih.gov/osra/

The wrapper comes with an automatic installation of all dependencies through the
galaxy toolshed.
