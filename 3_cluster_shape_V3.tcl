set file [read [open "v1_I_simple.dat" r]]
set output [open "v1_simple_I_shape.dat" w]

proc sum {args} {
  set result 0.0
  foreach n $args {
     set result [expr $result + $n]
  }
  return $result
}

set x {}
set y1 {}
set y2 {}
set y3 {}

puts $output "#SIZE NUMBEROFCLUSTERS AVERAGERATIO SUM STD "

for {set size 1} {$size <= 240} {set size [expr $size + 1 ]} {

set shapedat {}
### shapedat2 is used to calcuate STD
set shapedat2 {}     

foreach {a b c d} $file {
#puts $output2 "$a $b $c $d $e $f"
if { [expr $size] == $b} {

lappend shapedat $c
lappend shapedat2 [expr $c*$c]
#puts "FRAME $a CLUSTER $b "

}
}

# This if{} will add 0 to the array and increases the count by 1 unit. 

if { [llength $shapedat] == 0.0} {
lappend shapedat 0.0
#lappend shapedat2 0.0
}

puts "SIZE $size NUMBEROFCLUSTERS [expr [llength $shapedat]] AVERAGERATIO [expr [eval sum $shapedat]/[llength $shapedat]] SUM [eval sum $shapedat] STD [expr sqrt([eval sum $shapedat2]/[llength $shapedat2] -[eval sum $shapedat]*[eval sum $shapedat]/[llength $shapedat]/[llength $shapedat]) ] "

puts $output "$size [expr [llength $shapedat]] [expr [eval sum $shapedat]/[llength $shapedat]] [eval sum $shapedat] [expr sqrt([eval sum $shapedat2]/[llength $shapedat2] -[eval sum $shapedat]*[eval sum $shapedat]/[llength $shapedat]/[llength $shapedat]) ] "

#puts $output2 " $size [expr [llength $shapedat]] [expr sqrt([eval sum $shapedat2]/[llength $shapedat2] -[eval sum $shapedat]*[eval sum $shapedat]/[llength $shapedat]/[llength $shapedat]) ]"

#lappend x [expr $size ]
#lappend y1 [expr [llength $shapedat]]
#lappend y2 [expr [eval sum $shapedat]/[llength $shapedat]]
#lappend y3 [expr sqrt([eval sum $shapedat2]/[llength $shapedat2] -[eval sum $shapedat]*[eval sum $shapedat]/[llength $shapedat]/[llength $shapedat]) ] 

unset shapedat
#unset shapedat2
}

#puts "black line is Number of Clusters over Time, green line is Sum of ratio over Time"

#package require multiplot
#set plot2 [multiplot -x $x -y $y1 -xlabel Time(ns) -ylabel Number_of_Clusters -plot]
#$plot2 add $x $y2 -ylabel NumberOfClusters/Sumratio(AverageRatio) -linecolor green -plot
#$plot2 add $x $y3 -linecolor red -plot
