
Tools for analyzing macromolecular structure data, with emphasis on
protein - DNA interactions.

###############################################################################
makeProteinArray.py
###############################################################################

Preprocessing of single mmCIF file to extract protein chains and store coordinate
data in array format designed for easy use in structure analysis.

###############################################################################
makeDnaArray.py
###############################################################################

Preprocessing of singe mmCIF file to extract DNA chains and store coordinate
data in array format designed for easy use in structure analysis.

###############################################################################
contactMatrixByResType.py
##############################################################################

###############################################################################
distanceMatrixBySeq.py
##############################################################################

Calculates and plots as a heat map the distance matrix for two sets of
structures. input files must be .npz files, with keys that match those in 
inputs below.

coordinate arrays must have shape (N,maxResSize,3) where N is number of 
residues, maxResSize is number of atoms/res. Therefore, indices indicate 
(residue #, atom #, coordinate).
