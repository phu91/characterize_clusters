set file [read [open "v1_I_simple.dat" r]]
set output [open "v1_simple_I_count.dat" w]

proc sum {args} {
  set result 0
  foreach n $args {
     set result [expr {$result + $n}]
  }
  return $result
}

set x {}
set y1 {}
set y2 {}

#set AVGCLUSTERSIZE 0.0
puts $output "#FRAME NUMBEROFCLUSTERS SUM AVGCLUSTERSIZE"

for {set frame 1} {$frame <= 500 } {set frame [expr $frame + 1 ]} {

set clusterdat {}

foreach {a b r d} $file {

if {$frame == $a} {

lappend clusterdat $b

#puts "a=$a b=$b"
}
}
#puts "[eval sum $clusterdat]"

#puts "NUMBEROFCLUSTERS [expr [llength $clusterdat]]"


#puts $output "clusterdat $clusterdat"
puts "FRAME $frame NUMBEROFCLUSTERS [expr [llength $clusterdat]] SUM [eval sum $clusterdat]  AVGCLUSTERSIZE [expr {[eval sum $clusterdat] / [llength $clusterdat]}]"
puts $output "$frame [expr [llength $clusterdat]] [eval sum $clusterdat] [expr [eval sum $clusterdat] / [llength $clusterdat]]"


#lappend x [expr $frame * 0.95]
#lappend y1 [expr [llength $clusterdat]]
#lappend y2 [expr [eval sum $clusterdat] / [llength $clusterdat]]
#lappend y3 [eval sum $clusterdat]

unset clusterdat
}

#puts "black line is Number of Clusters over Time, blue line is Average Cluster Size over Time, NO red line is Sum of Cluster over Time"

#package require multiplot
#set plot2 [multiplot -x $x -y $y1 -xlabel Time(ns) -ylabel Number_of_Clusters -linecolor black -plot]
#$plot2 add $x $y2 -ylabel AverageClusterSize -linecolor blue -plot
#$plot2 add $x $y3 -ylabel SumClusterSize -fillcolor red -plot

#close $file
close $output
