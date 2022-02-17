# CW_domain_paper
Code to calculate the angle between two trp residues

## `01_define_input.R`
This script was used to collect all pdb structures and to define an `input_table.tsv` table that was used as an input parameter file for the `02_compute_trp_angles.R`. 


## `02_compute_trp_angles.R`
This script is used to calculate the angles between the truptophan rings.
A tryptophan ring was in this computation defined by the xyz coordinates of the atoms c("CG", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3", "CD2") from each trp residue. A 2d plane was fitted to the coordinates of the tryptophan ring with a linear model `lm(z ~ x + y)`.
The dihedral angle between two planes from two tryptophans was calculated from the angle between the normal vector to the two planes.
The angle and the 180-angle was computed and the correct orientation was selected manually from 3D plots
Only the angle from the A chain was used in the paper.
