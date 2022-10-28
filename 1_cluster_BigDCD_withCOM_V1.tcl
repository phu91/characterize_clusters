############### CLUSTER ALGORITHM ##################################
####################################################################
### Original source by Unknown
### Edited by Phu K. Tang
# Version 5.1
# 
# Notes:
# Set all to represent a whole mCPT #. Because mCPT has 12 residues
# thus to find the number of mCPT in a cluster, we need to divide
# $nagg by 12 to obtain the correct numbers of aggregations.
#
#	DATE			Notes
# 	Apr 13th, 2018		Add bigdcd.tcl
# 	Feb 22nd, 2019		Orient.tcl needs to return $s1, $s2, $s3



#####################################################################

############### REQUIREMENT #################
#############################################
source bigdcd.tcl			
package require Orient		
namespace import Orient::orient

############### INPUTS ######################
#############################################
#
#

#set mol [mol new "selection_beta.pdb" type pdb]
#mol addfile "sample.dcd" waitfor all 

#set resname1 "QTA"
#set drug1 "name C320 to C345 N333 N339 O324 and not type cA"
#set drug2 "name C424 to C451 N439 N445 O430 and not name C43 C44"
#set drug3 "name C370 to C398 N386 N392 O377 and not name C38 C39"
#set drug4 "name C265 to C292 N280 O271 N286"
#set all "$drug1 or $drug2 or $drug3 or $drug4"

#set popc [atomselect top "resname PA PC OL"]

### Pair-wise parameter ###
set rcut 3.0

### Get the residue number for each residue
set surf [atomselect top "type nb"]
set reslist [lsort -unique -real [$surf get beta]]
#lremove $reslist 0.0

$surf delete

### Get the numbers of frames. But we use BigDCD so this can be opted out.
#set nf [molinfo top get numframes]

### DCD file ###
set dcd "/home/phu/Desktop/tt/tt1_meta/meta-amber/72_random/eq-310K/noN-EQ-54ns-310K.dcd"
#set dcd "/home/phu/Desktop/da/qcpt/full_trajectory/II-short.dcd"
#set dcd "../25_us_qcpt_PBC_FIXED.dcd"
#set dcd "sample.dcd"

### Execution ###
bigdcd naggs $dcd
#bigdcd_wait

### Output
#set ofilena [open "COM_numbersOfAgg_25us_qCPT.dat" w]
#set ofilena [open "qcpt_cluster_full.dat" w]
set ofilena1 [open "v1_I_simple.dat" w]
set ofilena2 [open "v1_I_detail.dat" w]

############### PROCEDURE ###################
#############################################

proc naggs {frame} {
	global rcut reslist ofilena1 ofilena2
	
##########################
#### Fix PBC condition from Anton2 ####

set popc [atomselect top "type nb"]

pbc wrap -compound fragment -centersel $popc

#set all [atomselect top all]
#set gec [measure center $all]
#$all moveby [vecscale -1.0 $gec]
#$all delete  

#### Do BOUNDARY ATOMS ####

	### Create a smaller box inside the simulated box
    set dimx [expr [molinfo top get a frame $frame] / 2.0]
    set dimy [expr [molinfo top get b frame $frame] / 2.0]
    set dimz [expr [molinfo top get c frame $frame] / 2.0]
    set edgex [expr $dimx - $rcut]
    set edgey [expr $dimy - $rcut]
    set edgez [expr $dimz - $rcut]
    ### Get all the residues at the boundary
    set boundary [atomselect top "beta $reslist and not (z < $edgez and z > -$edgez and x < $edgex and x > -$edgex and y < $edgey and y > -$edgey)" frame $frame]

    set tmplist [lsort -unique -real [$boundary get beta]]
    $boundary delete
 
#     puts "TMPLIST ==> $tmplist  END"
foreach tmpres $tmplist {
	set tmp2 [atomselect top "beta $tmpres" frame $frame]
	set tmporig [measure center $tmp2]
	set xshift($tmpres) [lindex $tmporig 0]
	set yshift($tmpres) [lindex $tmporig 1]
	set zshift($tmpres) [lindex $tmporig 2]

	if { [lindex $tmporig 0] > [expr $dimx-$rcut] || [lindex $tmporig 0] < -[expr $dimx-$rcut] } {
	    if { $xshift($tmpres) < 0 } {
		set xshift($tmpres) [expr [lindex $tmporig 0] + [expr $dimx*2.0]]
	    } else {
		set xshift($tmpres) [expr [lindex $tmporig 0] - [expr $dimx*2.0]]
	    }
	}

	if { [lindex $tmporig 1] > [expr $dimy-$rcut] || [lindex $tmporig 1] < -[expr $dimy-$rcut] } {
	    if { $yshift($tmpres) < 0 } {
		set yshift($tmpres) [expr [lindex $tmporig 1] + [expr $dimy*2.0]]
	    } else {
		set yshift($tmpres) [expr [lindex $tmporig 1] - [expr $dimy*2.0]]
	    }
	}

	if { [lindex $tmporig 2] > [expr $dimz-$rcut] || [lindex $tmporig 2] < -[expr $dimz-$rcut] } {
	    if { $zshift($tmpres) < 0 } {
		set zshift($tmpres) [expr [lindex $tmporig 2] + [expr $dimz*2.0]]
	    } else {
		set zshift($tmpres) [expr [lindex $tmporig 2] - [expr $dimz*2.0]]
	    }
	}

	$tmp2 delete
    }
##### Find numbers of clusters #####

    set ncluster 0

    # Notes: set reslist [lsort -unique -real [$surf get residue]] (--above--)
    #### This loop finds the neighbors in a cluster ####
    foreach res $reslist {
    	### Set neighbors to be anyone who is within the rcut distance
		set neighbor [atomselect top "beta $reslist and within $rcut of beta $res" frame $frame]
		set nlist [lsort -real -unique [$neighbor get beta]]
        $neighbor delete
		set cluster($ncluster) $nlist
		set tmp1 [atomselect top "beta $res" frame $frame]
	#### Compare between the neighbor list and all selected residues ####
	foreach tmpres $tmplist {
	    if { [veclength [vecsub [list $xshift($tmpres) $yshift($tmpres) $zshift($tmpres)] [measure center $tmp1]]] < $rcut } {
		lappend cluster($ncluster) $tmpres
#		puts $cluster($ncluster)
	    }
	}
	$tmp1 delete
	incr ncluster
    }

    puts "CLUSTERS: $ncluster"
    		######		THE END OF CLUSTER ALGORITHM	#####

############### ORIENTATION #################
#############################################
	#	Loop over each cluster and align it to the principal axises
    #### while {0} = stop the loop; {1} = keep running the loop ####
    set loop 1
    while { $loop } {
	set loop 0
#	puts "COMBINE $ncluster"
	for {set i 0} {$i < [expr $ncluster - 1]} {incr i} {
	    for {set j [expr $i +1]} {$j < $ncluster} {incr j} {
		set found 0
		foreach ni $cluster($i) {
		    foreach nj $cluster($j) {
			if { $ni == $nj } {
			    set found 1
                	    set loop 1   
			}
		    }
		}
		if { $found == 1 } {
		    set cluster($i) [lsort -unique -real [concat $cluster($i) $cluster($j)]]
		    for {set k $j} {$k < [expr $ncluster - 1] } {incr k} {
			set cluster($k) $cluster([expr $k + 1])
		    }
		    set j [expr $j - 1] 
		    set ncluster [expr $ncluster - 1]
		    #		    puts "NCLUSTERS: $ncluster"
		}
		}
	}
#	puts $ncluster
    }

    puts "FOUND $ncluster CLUSTERS"
#    set ccount 0
    set ntot 0
    for {set i 0} {$i < $ncluster} {incr i} {

	set nagg [llength $cluster($i)]
	set ntot [expr $ntot + $nagg]
	if {$nagg > 1} {


		foreach ni $cluster($i) {
		set member $ni
	#	puts $OUT2 "$ni"
		set tmp2 [atomselect top "beta $ni" frame $frame]
		set tmporig [measure center $tmp2]
	if { [lindex $tmporig 0] > [expr $dimx-$rcut] || [lindex $tmporig 0] < -[expr $dimx-$rcut] } {
	    if { $xshift($tmpres) < 0 } {
		set xshift($tmpres) [expr [lindex $tmporig 0] + [expr $dimx*2.0]]
	    } else {
		set xshift($tmpres) [expr [lindex $tmporig 0] - [expr $dimx*2.0]]
	    }
	}

	if { [lindex $tmporig 1] > [expr $dimy-$rcut] || [lindex $tmporig 1] < -[expr $dimy-$rcut] } {
	    if { $yshift($tmpres) < 0 } {
		set yshift($tmpres) [expr [lindex $tmporig 1] + [expr $dimy*2.0]]
	    } else {
		set yshift($tmpres) [expr [lindex $tmporig 1] - [expr $dimy*2.0]]
	    }
	}

	if { [lindex $tmporig 2] > [expr $dimz-$rcut] || [lindex $tmporig 2] < -[expr $dimz-$rcut] } {
	    if { $zshift($tmpres) < 0 } {
		set zshift($tmpres) [expr [lindex $tmporig 2] + [expr $dimz*2.0]]
	    } else {
		set zshift($tmpres) [expr [lindex $tmporig 2] - [expr $dimz*2.0]]
	    }
	}
		set xshift($tmpres) [lindex $tmporig 0]
		set yshift($tmpres) [lindex $tmporig 1]
		set zshift($tmpres) [lindex $tmporig 2]

		}
		#    orient the cluster along xy plane, write out to separate file
		set sel [atomselect top "beta $cluster($i)" frame $frame]
		set countlist [lsort -real -unique [$sel get beta]]		
		puts $ofilena2 "MEMBERS IN CLUSTER===> $countlist"
### Export residues in Templist

		set commb [lindex [measure center $popc] 2]
                set pos2 [lindex [measure minmax [atomselect top "beta $countlist"]] 0 2]
                set dz [vecdist $pos2 $commb]
		
		set pos [measure center $sel]
		set xpos [lindex $pos 0]
		set ypos [lindex $pos 1]
		set zpos [lindex $pos 2]
		set xtran [expr -1.0*$xpos]
		set ytran [expr -1.0*$ypos]
		set ztran [expr -1.0*$zpos]
		$sel moveby "{$xtran} {$ytran} {$ztran}"

### Calculate the moment of Inertia (I) using Orient package. CHECK THE Orient.tcl
#	for more information. $s1, $s2, $s3 from the Orient.tcl 
#	will hold 3 components of the diagonalized I. They must be returned from
#	the Orient.tcl. The new Orient.tcl may not have it. 
#
#	$maj1 = Imax and $maj3 = Imin.  
#
		set I [draw principalaxes $sel]
  		set A [orient $sel [lindex $I 2] {0 0 1}]
  		$sel move $A
		set I [draw principalaxes $sel]
  		set A [orient $sel [lindex $I 1] {0 1 0}]
		set maj1 [lindex $I 3]
		set maj2 [lindex $I 4]
		set maj3 [lindex $I 5]

#		foreach ni $cluster($i) {
#		set j [expr $j + 1.0]
#		set member $ni
#	#	puts $OUT2 "$ni"
#		set tmp2 [atomselect top "residue $ni" frame $frame]
#		set type [$tmp2 get resname]
#		set type1 [lindex $type 0]
#		}

		set shape [expr $maj3/$maj1]
  		$sel move $A

### Center of mass

#		set commb [lindex [measure center $popc]]
#		set pos2 [lindex [measure center [atomselect top "beta $countlist"]]]
#		set dz [vecdist $pos2 $commb]
### Printing Output ###
		
		puts $ofilena1 "$frame [expr [llength $countlist]] $shape $dz\n"
		puts $ofilena2 "FRAME = $frame SIZE = [expr [llength $countlist]] SHAPE = $shape MB = $commb CLUSTER = $pos2 DZ = $dz\n"
#	       	puts $ofilena1 "frame = $frame mb = $commb cluster = $pos2 dz = $dz\n"

### Cleaning the buffer to speed up the process.
	    	flush $ofilena1
		flush $ofilena2
		}
flush stdout

	}
}

