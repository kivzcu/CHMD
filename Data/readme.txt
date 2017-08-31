Structure of data
-----------------
REQUIRED:
m0.vtk		= musculotendon unit to decompose

f0.vtk
f1.vtk
..
fn.vtk		= superficial musculotendon fibres lying on the surface of the muscle, each polyline represents one fibre, could be in one file or multiple files,
		  if points in a file are unconnected, they are connected automatically in one fibre

ab0.vtk
ab1.vtk
..
abn.vtk		= bones with origin attachment areas of the mtu

bb0.vtk
bb1.vtk
..
bbn.vtk		= bones with insertion attachment areas of the mtu

OPTIONAL:
ao0.vtk
ao1.vtk
..
aon.vtk		= surface mesh of the origin attachment area

bo0.vtk
bo1.vtk
..
bon.vtk		= surface mesh of the insertion attachment area


	