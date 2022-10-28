# characterize_clusters
These are TCL scripts that are used in tk Console of VMD. **Orient package is required and may not be supported in the future VMD.** 

1. `source bigdcd.tcl`                        : call BIGDCD procedure **There are room to optimize the script for better memory usgae and speed**
2. `source  1_cluster_BigDCD_withCOM_V1.tcl`  : Defining cluster. **BUGS potential for multi-resname monomer !! Use with cautions!!**
3. `source 2_cluster_count_V3.tcl`            : COUNTING found clusters
4. `source 3_cluster_shape_V3.tcl`            : CHARACTERIZING THE SHAPE of found clusters
