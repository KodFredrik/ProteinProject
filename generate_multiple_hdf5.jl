# EXAMPLE 1
import PyCall
import PyPlot
PyCall.@pyimport matplotlib
import FFTW
import cidsim
#Only generates first file in folderIn
function generate_multiple(folderIn, folderOut)
	@time begin
	# INPUT PARAMETERS
	BEAM_ENERGY = 8265.647
	PULSE_ENERGY = 1E-3 # J
	FOCUS_DIAMETER = 1E-6 # m
	BANDWIDTH = 0.0
	DETECTOR_DISTANCE = 0.1
	PIXEL_SIZE = 500E-6
	DETECTOR_ARRAY = (550, 550)
	POLARIZATION = (1, 0)
	BANDWIDTH_SLICES = 1
	QUANTUM_EFFICIENCY = 1.0
	NUFFTTOL = 1e-10
	NTHREADS = Threads.nthreads()
	    
	#EDIT PDB FILE LOCATION HERE
	foldern = readdir(folderIn)
	for filename in foldern
		PDB_FILE = folderIn*"/"*filename
		# DO DIFFRACTION
		s = cidsim.Source.flat(BEAM_ENERGY, PULSE_ENERGY, FOCUS_DIAMETER, BANDWIDTH, 			BANDWIDTH_SLICES, POLARIZATION)
		d = cidsim.Detectors.detector(DETECTOR_DISTANCE, PIXEL_SIZE, DETECTOR_ARRAY, QUANTUM_EFFICIENCY)
		pdb = cidsim.Target.pdb(PDB_FILE)
		result = cidsim.Diffraction.scattering(s, d, pdb, "coherent", "tabulated", NUFFTTOL, NTHREADS)
		# SAVE TO FILE, EDIT OUTPUT HDF5 LOCATION
		cidsim.write2hdf5(folderOut*"/"*filename[1:end-4]*".h5", s, d, result)
	end
	end
end

infolder = "removed_noise_pdb"
outfolder = "removed_hdf5"
generate_multiple(infolder,outfolder)

#Translations
#trs_in = ["pdb_441_trsl_0.5Ang_noise", "pdb_441_trsl_1Ang_noise", 
#"pdb_441_trsl_5Ang_noise", "pdb_441_trsl_10Ang_noise"]
#trs_out = ["trsl_0.5Ang_noise", "trsl_1Ang_noise", "trsl_5Ang_noise", "trsl_10Ang_noise"]
#for i in 1:length(trs_in)
#	generate_single(trs_in[i],trs_out[i])
#end

#Rotations
#rot_in = ["pdb_data_rot1_noise", "pdb_data_rot5_noise", "pdb_data_rot10_noise", "pdb_data_rot20_noise", "pdb_data_rot45_noise",]
#rot_out = ["rot1_noise", "rot5_noise", "rot10_noise", "rot20_noise", "rot45_noise"]
#for i in 1:length(rot_in)
#	generate_single(rot_in[i],rot_out[i])
#end
