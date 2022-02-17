# This script was used to download and prepare for angle calculations


# Using 
# R 4.1.2 for windows 64 bit


# rstudioapi
# bio3d
# httr




library("bio3d")

# Define working dir relative to this document
cw = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(cw)


# Read input table
# Contains proteins that should be analyzed
input_table = read.delim(file = "input_table.tsv")

# Download all protein structures
# Only do once
# dir.create("pdb_files",showWarnings = F)
# setwd("./pdb_files")
# sapply(input_table$PDB_ID, get.pdb)
# setwd(cw)

# print
#2L7P         6QXZ         5YVX         4QQ4         5SVY         5SVX         5IX2         5IX1         2E61         2RR4         4O62 
#"./2L7P.pdb" "./6QXZ.pdb" "./5YVX.pdb" "./4QQ4.pdb" "./5SVY.pdb" "./5SVX.pdb" "./5IX2.pdb" "./5IX1.pdb" "./2E61.pdb" "./2RR4.pdb" "./4O62.pdb" 

# Find the location of all trp residues in the CW domains
x = input_table$PDB_ID[1]




#input_row = input_table[4,]
trp_distances = apply(input_table, 1, function(input_row){ # loop over each 
  pdb_id = input_row[2]
  pdb_chains = strsplit(as.character(input_row[3]), split = ",")[[1]]
  
  pdb = read.pdb(paste0("./pdb_files/", pdb_id , ".pdb"))
  
  # Find tryptophan residues and their CA atom
  coords = pdb$atom[pdb$atom$resid == "TRP",]
  
  # Keep only relevant chains
  coords = coords[coords$chain %in% pdb_chains,]
  
  # Keep only relevant atoms
  coords = coords[coords$elety %in% "CA",]
  rownames(coords) <-paste(coords$chain, coords$resno)
  
  
  # Split by chain
  coords = split(coords, f = coords$chain)
  
  lapply(coords, function(chain){
    dist(chain[c("x", "y", "z")])
  })
})

names(trp_distances) <- input_table$PDB_ID

trp_distances
  
  

  