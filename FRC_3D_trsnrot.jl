import PyCall
import PyPlot
PyCall.@pyimport matplotlib
import FFTW
import cidsim
import HDF5
using Plots; pyplot()

function mu_pixel_avg(foldername)
	path = "/home/fredrik/Skrivbord/PROJDIR/"*foldername*"/"
	set_folder = readdir(foldername)
	#Initiate average intensity matrix mu1
	global mu_sum = zeros(Float64,550,550)
	for filename in set_folder
		local fid = HDF5.h5open(path*filename, "r")
		#Find dataset within groups
		local entry_1 = fid["entry_1"]
		local data_1 = entry_1["data_1"]
		local data = data_1["data"]
		#Read dataset
		local intensity = read(data)
		local mu_new = mu_sum .+ intensity
		global mu_sum = mu_new
		close(fid)
	end
	#Divide mu1_sum by number of matrices in set1
	M = length(set_folder)
	mu = mu_sum./M
	return mu
end

function hdf5_to_data(foldername,filename, dataset)
	path = "/home/fredrik/Skrivbord/PROJDIR/"*foldername*"/"
	set_folder = readdir(foldername)
	fid = HDF5.h5open(path*filename, "r")
	entry_1 = fid["entry_1"]
	data_1 = entry_1["data_1"]
	data = read(data_1[dataset])
	close(fid)
	return data
end

function mu_rad_profiles(mu,smatrix)
	#Define radial profile (Sebastians profile below)
	RADIALSTEP = smatrix[1,1]/200
	sr = collect(0:RADIALSTEP:smatrix[1,1])
	profiles_list = Array{Any}(undef, length(sr)-1)
	# final criteria for selection
	for i in 1:length(sr)-1
		selectedarea = (sr[i].<smatrix) .* (smatrix.<=sr[i+1])
		profiles_list[i] = selectedarea
	end
	SR = Vector{Float64}(undef,200)
	for j in 1:length(sr)-1
		SR[j] = (sr[j]+sr[j+1])/2
	end
	return profiles_list, SR
	end

function mu_radial_average(mu, profiles_list)
	mu_ravg = zeros(0)
	for i in 1:length(profiles_list)
		append!(mu_ravg,sum(mu[profiles_list[i]])/length(mu[profiles_list[i]]))
	end
	return mu_ravg
end

##should get a radial profile to pass in
function FRC(profiles ,mu1 ,mu2,mu1_radavg,mu2_radavg, SR)
	frc_list = zeros(0)
	for i in 1:length(profiles)
		factor1 = sum((mu1[profiles[i]] .- mu1_radavg[i]) .* (mu2[profiles[i]] .- mu2_radavg[i]))
		factor2 = sqrt(sum((mu1[profiles[i]] .- mu1_radavg[i]).^2))
		factor3 = sqrt(sum((mu2[profiles[i]] .- mu2_radavg[i]).^2))
		frc = factor1 / (factor2*factor3)
		if isnan(frc) == true
			frc = 0
		end
		append!(frc_list, frc)
	end
	return frc_list
end

function interpolate(x,y,z)
	lx = length(x)
	ly = length(y)
	init = zeros(ly,ly)
	int1rv = x*(length(y)*(1/x[end]))
	xnew = collect(LinRange(0,x[end],length(y)+1))
	xnew = xnew[2:end]
	plane = ones(200,200)/2
	#plane[:,1] = z[:,1]
	init[:,1] = z[:,1]
	for i in 1:(length(x)-1)
		interv = length(y)*(x[i]/x[end])
		k = (z[:,i+1].-z[:,i])/(int1rv[i+1]-int1rv[i])
		n = 0
		if i==1
			for j in 1:floor(Int,int1rv[i+1])
			init[:,j] = z[:,i] + n*k
			n += 1
			end
			else
				for j in floor(Int,int1rv[i]):floor(Int,int1rv[i+1])
				init[:,j] = z[:,i] .+ n*k
				n += 1
				end
		end
		#plane[:,floor(Int,int1rv[i])] = z[:,i]
	end
	init[:,floor(Int,int1rv[end])] = z[:,end]
	#plane[:,floor(Int,int1rv[end])] = z[:,end]
	return xnew ,y , init, plane
end

###########INPUT###########
#Titel = "Translation and 0.2 Å std noise" #TRS TITEL
Titel = "Rotation and 0.2 Å std noise" #ROT TITEl
#Titel = "Translation and 10*rotation and 0.2 Å std noise with host"
savetitel = "ROT_FINAL.png"
#Xlabel = "Translation [Å]"
Xlabel = "Degrees"
#Xlabel = "Translation and 10*rotation"
infolder = "ROT_FINAL_HDF5"

sx = hdf5_to_data(infolder, "0_0.h5", "sx")
sy = hdf5_to_data(infolder, "0_0.h5", "sy")
s = 2 * pi * 1E-10 * sqrt.(sx.^2 + sy.^2)

###Noise 3D plots
#noise = [1, 5, 10, 15, 20, 25, 50, 100]
#trsstr = ["0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1", "1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "1.9", "2", "2.1", "2.2", "2.3", "2.4", "2.5", "2.6", "2.7" ,"2.8", "2.9", "3" ,"3.1", "3.2", "3.3", "3.4", "3.5", "3.6", "3.7", "3.8", "3.9", "4", "4.1", "4.2", "4.3", "4.4", "4.5", "4.6", "4.7", "4.8", "4.9", "5"]
#rotstr = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"]
#rot = [0, 1, 2, 3 ,4 ,5, 6, 7, 8, 9, 10, 11, 12, 13 ,14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]
#trs = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5 ,1.6 ,1.7 ,1.8 ,1.9 ,2 ,2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9 ,5]
rotrange = LinRange(0,60,60)
floatrange = [round(i) for i in rotrange]
rot = Int.(floatrange)
rotstr = [string(i) for i in rot]
trsrange = LinRange(0,8,80)
trs = [round(i, digits=1) for i in trsrange]
trsstr = [string(i) for i in trs]
#Determine what to compare with, here noised ideal with std ~0.1 Å
#EDIT here and also output of interpolation
mu1 = hdf5_to_data(infolder*"/","0_0.h5","data")
#combo = [trsstr[1]*"_"*j*".h5" for j in rotstr] #use for rotation
combo = ["0"*"_"*j*".h5" for j in rotstr] #use for emergencies
#combo = [i*"_"*rotstr[1]*".h5" for i in trsstr] #use for translation
#combo = [trsstr[i]*"_"*rotstr[i]*".h5" for i in 1:length(trsstr)] #use for diagonal
#Titel = "Translation and 0.2 Å std noise"
#savetitel = "HOST.png"
#Xlabel = "Translation and 10*Degrees"

#PyPlot.figure()
FRC_LIST = zeros(200,length(combo))
R_LIST = zeros(length(combo))
for i in 1:length(combo)
	mu2 = hdf5_to_data(infolder, combo[i],"data")
	global profiles, SR = mu_rad_profiles(mu2, s)
	mu1_radavg = mu_radial_average(mu1, profiles)
	mu2_radavg = mu_radial_average(mu2, profiles)
	FRC_list = FRC(profiles ,mu1 ,mu2,mu1_radavg,mu2_radavg, SR)
	FRC_LIST[:,i] = FRC_list
	#Implement R-factor
	local Iideal = mu1
	local Ireal = mu2
	local K = sum(sqrt.(Ireal))/sum(sqrt.(Iideal))
	local R = abs((1/K)*sum(sqrt.(Ireal)-sqrt.(Iideal))/sum(sqrt.(Iideal)))
	local RoundR = string(round.(R; digits=2))
	R_LIST[i] = round.(R; digits=2)
	#Edit this for every change
	#PyPlot.plot(SR, FRC_list, label=rotstr[i]*"Deg")
	end

#R-factor plot
PyPlot.figure()
PyPlot.title("R-factor")
#PyPlot.xlabel("Translation [Å]")
PyPlot.xlabel("Degrees")
PyPlot.ylabel("R-factor")
PyPlot.plot(rot, R_LIST)
PyPlot.savefig("R-factor")

#####2DPlot######
#PyPlot.title(Titel)
#PyPlot.xlabel("q [1/Å]")
#PyPlot.ylabel("FRC")
#PyPlot.legend(loc="lower left")
#PyPlot.savefig("2D"*Titel*".png")

##get radial profile plot
#profiles, SR = mu_rad_profiles(mu1, s)
#mu1_radavg = mu_radial_average(mu1, profiles)
#
#PyPlot.figure()
#PyPlot.imshow(profiles[100])
#PyPlot.title("Radial profile (200 steps)")
#PyPlot.savefig("Radial_profile_200steps.png")
#Plot ideal (mu1) radial average
#PyPlot.figure()
#PyPlot.plot(SR, mu1_radavg)
#PyPlot.title("Radial average of a 4x4x1 1mr0 crystal with noise")
#PyPlot.xlabel("q [1/Å]")
#PyPlot.yscale("log")
#PyPlot.ylabel("Radial Average")
#PyPlot.savefig("ravg.png")



#plane = ones(200,200)/2
#EDIT
x,y,z,plane = interpolate(rot,SR,FRC_LIST)

PyPlot.figure(figsize=(10,10))
s = PyPlot.plot_surface(x, y, z, cmap="seismic",edgecolor="none",linewidth=1,alpha=0.8,vmin=0,vmax=1)
PyPlot.plot_surface(x, y, plane, cmap="viridis",edgecolor="none",linewidth=1,alpha=0.8)
PyPlot.title(Titel)
PyPlot.ylabel("q [1/Å]")
PyPlot.xlabel(Xlabel)
cbar1 = PyPlot.colorbar(s, pad=0.1)
cbar1.set_label("FRC")
PyPlot.view_init(5,45,0)
#EDIT
PyPlot.savefig("3D"*savetitel)

#overindex = findall(x->x > 0.5, z)
#PyPlot.figure()
#PyPlot.imshow(z, cmap="seismic", interpolation= "nearest", extent = [x[1],x[end],y[end],y[1]], aspect="auto" )
#PyPlot.colorbar()
#PyPlot.ylabel("q [1/Å]")
#PyPlot.xlabel(Xlabel)
#PyPlot.title(Titel)
#PyPlot.savefig("heat"*savetitel)

#PyPlot.close()
PyPlot.figure()
PyPlot.imshow(z', cmap="seismic", interpolation= "nearest", extent = [y[1],y[end],x[end],x[1]], aspect="auto" ,vmin=0,vmax=1)
cbar = PyPlot.colorbar()
cbar.set_label("FRC")
PyPlot.ylabel(Xlabel)
PyPlot.xlabel("q [1/Å]")
#PyPlot.zlabel("FRC")
PyPlot.title(Titel)
PyPlot.savefig("heat_"*savetitel)
PyPlot.show()
