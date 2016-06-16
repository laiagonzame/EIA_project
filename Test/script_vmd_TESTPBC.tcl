mol delete top
mol load xyz test_pbc.xyz
mol delrep 0 top
dislay resetview

set selC[atomselect top ""name C""]
$selC set raidus 0.01
color Display Background white
mol representation VDW 0.20000 16.0000
mol selection name C
mol material Opaque
mol color ColorID 1
mol addrep top

pbc set {4 4 4} -all -molid top
pbc box -center origin
