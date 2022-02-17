# compute angles between trp ring planes


# Using
# R 4.1.2 for windows 64 bit
# rstudioapi
# bio3d
# hfit
# hyper.fit

library("hfit")
library("hyper.fit")
library("bio3d")

# Define working dir relative to this document
cw = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(cw)

# Define a trp plane as the atoms in the trp ring 
def_trp_ring_atoms = c("CG", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3", "CD2")

# Functions
# This function filters a dataframe and returns all the atoms in a pdb file that matches the creteria
getAtoms = function(pdb_id, pdb_chains, pdb_resno, def_trp_ring_atoms){
  pdb = read.pdb(paste0("./pdb_files/", pdb_id , ".pdb"))
  # Keep only relevant residue numbers
  coords = pdb$atom[pdb$atom$resno %in% pdb_resno,]
  
  # Keep only relevant chains
  coords = coords[coords$elety %in% def_trp_ring_atoms,]
  
  #coords_split = split(coords, f = coords$resno)
  return(coords)
}
# test getCoords
# single
# i = 4
# getAtoms(pdb_id = input_table$PDB_ID[i], pdb_chains = input_table$Chain[i] , pdb_resno = strsplit(input_table$Trp_number[i], split = ",")[[1]], def_trp_ring_atoms = def_trp_ring_atoms )
# Whole data.frame
# lapply(1:nrow(input_table), function(i){
#   getAtoms(pdb_id = input_table$PDB_ID[i], pdb_chains = input_table$Chain[i] , pdb_resno = strsplit(input_table$Trp_number[i], split = ",")[[1]], trp_ring_atoms = trp_ring_atoms )
#   
# })



# Computes angles between planes
# input should be a dataframe with a, b, c, and d columns
getAngleBetweenPlanes = function(x){
  cos_theta = (x$a[1]*x$a[2] + x$b[1]*x$b[2] + x$c[1]*x$c[2])/(sqrt(x$a[1]^2 + x$b[1]^2 + x$c[1]^2) * sqrt(x$a[2]^2 + x$b[2]^2 + x$c[2]^2))
  angle1 = 180*acos(cos_theta)/pi
  angle2 = 180-angle1
  c(angle1, angle2)
}
# getAngleBetweenPlanes(chain)

# Run

# Read input table
# Contains proteins that should be analyzed
input_table = read.delim(file = "input_table.tsv")

# Extract atoms that define the trp ring 
i = 4
trp_ring_atoms = lapply(1:nrow(input_table), function(i){
  # Get atoms
  atoms = getAtoms(pdb_id = input_table$PDB_ID[i], pdb_chains = strsplit(input_table$Chain[i], split = ",")[[1]] , pdb_resno = strsplit(input_table$Trp_number[i], split = ",")[[1]], def_trp_ring_atoms = def_trp_ring_atoms )
  
  # split into pdb chains
  atoms = split(atoms, f = atoms$chain)
  atoms
})

# Fit a linear model to each trp ring
x = trp_ring_atoms[[4]]
trp_ring_planes = lapply(trp_ring_atoms, function(pdb){ # loop over pdb
  chain = pdb[[1]]
  lapply(pdb, function(chain){ # loop over chain 
    do.call(rbind, lapply(split(chain, f = chain$resno), function(coord_df){ # loop over residue number
      # Fit both rings to a linear model
      fit <- lm(coord_df$z ~ coord_df$x + coord_df$y)
      
      # Extract coefficients from fit
      coefs <- coef(fit)
      a <- coefs["coord_df$x"]
      b <- coefs["coord_df$y"]
      c <- -1
      d <- coefs["(Intercept)"]
      
      # Return data.frame with coefficients
      df = data.frame(a = a, b = b, c = c, d = d)
      row.names(df) = NULL
      df
    }))
  })
})

# Plotting and visual QC
i = 4 # pdb to plot
chain = "A" # Chain to plot
plot3d(trp_ring_atoms[[i]][[chain]]$x, trp_ring_atoms[[i]][[chain]]$y, trp_ring_atoms[[i]][[chain]]$z, type = "s", col = "red", size = 1)
planes3d(trp_ring_planes[[i]][[chain]]$a, trp_ring_planes[[i]][[chain]]$b, trp_ring_planes[[i]][[chain]]$c, trp_ring_planes[[i]][[chain]]$d, type = "s", col = "red", size = 1)


# compute angle betwee planes

pdb = trp_ring_planes[[7]]
trp_angles = lapply(trp_ring_planes, function(pdb){
  # chain = pdb[[1]]
  df = do.call(rbind,lapply(pdb, function(chain){
    sort(getAngleBetweenPlanes(chain))
  }))
  
  colnames(df) <- c("angle1", "angle2")
  df = as.data.frame(df)
  df$chain = rownames(df)
  df
})
names(trp_angles) = input_table$PDB_ID

# make table 


# Depending on the rotation of the molecule the "outer" or "inner" angle are both shown here.
# The angle corresponding closest to the manually measured one was selected.
trp_angles_df = do.call(rbind, trp_angles)
trp_angles_df

# manually selecting the inner angle
angle_selection_manual = c(1,2,2,2,2,2,2,1,1,1,1,2,2,2,2,2)
# set to zero if not selected

trp_angles_df$angle1[angle_selection_manual == 2] = 0
trp_angles_df$angle2[angle_selection_manual == 1] = 0
trp_angles_df$selected = trp_angles_df$angle1 + trp_angles_df$angle2



trp_angles_df



