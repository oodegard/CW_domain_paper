# compute angles between trp ring planes

library("bio3d")
library("hyper.fit")


# Define working dir relative to this document
cw = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(cw)

# Define a trp plane as the atoms in the trp ring 
trp_ring_atoms = c("CG", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3", "CD2")

# Read input table
# Contains proteins that should be analyzed
input_table = read.delim(file = "input_table.tsv")



# Test code
input_row = input_table[1,]
pdb_id = input_row[2]
pdb_chains = strsplit(as.character(input_row[3]), split = ",")[[1]]
pdb_resno = strsplit(as.character(input_row[4]), split = ",")[[1]]
pdb = read.pdb(paste0("./pdb_files/", pdb_id , ".pdb"))

# Keep only relevant residue numbers
coords = pdb$atom[pdb$atom$resno %in% pdb_resno,]
  
# Keep only relevant chains
coords = coords[coords$elety %in% trp_ring_atoms,]

coords_split = split(coords, f = coords$resno)


hyper.fit( coords_split[[1]][c("x", "y", "z")])

summary(hfit)

plot(hfit)
