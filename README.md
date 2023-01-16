# ProteinProject
The script translate_and_rotate.tcl is executed in VMD to generate the translated and rotated crystals.
The pdb files are moved into a folder in the Julia project directory.
The script generate_multiple_hdf5.jl is executed to simulate scattering and to generate the hdf5 files from the folder with the pdb
files to a new folder.
Execute the script FRC_3D_trsnrot.jl to generate the FRC plots. 
In this script mu1 is the reference image and mu2 is the "real" image. 
The script is not very nice to read but you can ask me at fredd3_95@hotmail.com is by any chance you would have to run it.
