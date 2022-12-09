#By the central limit theorem the sum of many random 
#variables will have an approximate normal distribution
#This procedure generates a random number between -0.5 and 0.5 from a normal distribution

proc norm_number {} {
	set x 0
	set y 10000
	set sigma [expr 1.0 / 12]
	set sigmaX [expr sqrt([expr $sigma / $y])]
	#Sum y random numbers between 0 and 1
	for {set i 0} {$i < $y} {incr i} {
		set x [expr $x + rand()]
	}
	#Normalize sum
	set transl [expr [expr $x/$y -0.5] / $sigmaX]
	return $transl
}
#
#NOTICE N HAS TO BE ALTERED DEPENDING ON CRYSTAL SIZE
#The input "scale" scales the standard deviation sigmaX by this factor in theory, which in this case is the sample variance for N samples.

proc translate {scale} {
	set N 16
	set pi 3.1415926535897931
	set movelist {}
	set variance 0
	for {set i 0} {$i < $N} {incr i} {
		#Select protein ID
		set fragID [expr $i]
		#Generate random number, scale, store in list for documentation
		set amount [norm_number]
		set scaledamount [expr $amount*$scale]
		lappend movelist [format "%.2f" $scaledamount]
		#Select a protein
		set sel [atomselect top "fragment $fragID"]
		#If we want to display rotation values in blue color on the proteins 8)
		#color Display Background white
		#set com [measure center $sel weight mass]
		#set dispamount [format %.2f $scaledamount]
		#graphics top text $com "$dispamount"
		#Randomize direction to translate in, then move molecule
		set number [expr rand()]
		set move_x [expr $scaledamount*cos(2*$pi*$number)]
		set move_y [expr $scaledamount*sin(2*$pi*$number)]
		$sel moveby [list $move_x $move_y 0]
		#Find variance 
		set variance [expr $variance + [expr $amount**2]]
	}
	set Variance [format "%.2f" [expr $variance/$N]]
	set std [format "%.2f" [expr sqrt($variance/$N)]]
	set scaled_std [format "%.2f" [expr $scale*sqrt($Variance)]]
	puts "The sample variance is $Variance times (SigmaX)^2"
	puts "The sample std is $std Å times SigmaX"
	puts "The scaled sample std is $scaled_std Å"
#	graphics top text {0 230 50} "Std = $scaled_std Angstrom"
	#return $movelist
#	return $Variance
}

proc rotate {stddegrees} {
	#N is the number of proteins
	set N 16
	set rotationlist {}
	set variance 0
	for {set i 0} {$i < $N} {incr i} {
		set fragID [expr $i]
		set sel [atomselect top "fragment $fragID"]
		#Randomize orientation, store angles, then rotate.
		set number [norm_number]
		lappend rotationlist [format "%.2f" [expr $number*$stddegrees]]
		set com [measure center $sel weight mass]
		#If we want to display rotation values in blue color on the proteins 8)
		set degnumber [format %.2f [expr $number*$stddegrees]]
		graphics top text $com "$degnumber"
		set matrix [transaxis z  [expr $number*$stddegrees]]
		$sel moveby [vecscale -1.0 $com]
		$sel move $matrix
		$sel moveby  $com
		#Find variance 
		set variance [expr $variance + [expr $number**2]]
	}
	set Variance [format "%.2f" [expr $variance/$N]]
	set std [format "%.2f" [expr sqrt($variance/$N)]]
	set scaled_std [format "%.2f" [expr $stddegrees*sqrt($Variance)]]
	puts "The sample variance is $Variance times (SigmaX)^2"
	puts "The sample std is $std degrees times SigmaX"
	puts "The scaled sample std is $scaled_std degrees"
#	return $Variance
}

set trs {3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5}
#{0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3}
#{0 0.1 0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.9 2.1 2.3 2.5 2.7 2.9 3.1 3.3 3.5 3.7 3.9 4.1 4.3 4.5 4.7 4.9}
# 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 4.4 4.6 4.8 5}
set rot {31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50}
#{0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30}
#1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45}
# 1 5 10 15 20 25 30 35 40 45}

#Create pdb files 
#for {set g 0} {$g < [llength $trs]} {incr g} {
#	for {set k 0} {$k < [llength $rot]} {incr k} {
#		set mol [mol load pdb 441_1mr0.pdb]
#		translate [lindex $trs $g]
#		rotate [lindex $rot $k]
#		animate write pdb [lindex $trs $g]_[lindex $rot $k].pdb $mol
#		mol delete all
#	}
#}

#Create diagonal elements
for {set g 0} {$g < [llength $trs]} {incr g} {
	set mol [mol load pdb 441_1mr0.pdb]
		translate [lindex $trs $g]
		rotate [lindex $rot $g]
		animate write pdb [lindex $trs $g]_[lindex $rot $g].pdb $mol
		mol delete all
}

