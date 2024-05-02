from solid2 import *

set_global_fn(100)

d = cube(5) + sphere(5).right(5) - cylinder(r=2, h=6)

# polyhedron(
#   points=[ [10,10,0],[10,-10,0],[-10,-10,0],[-10,10,0], // the four points at base
#            [0,0,10]  ],                                 // the apex point 
#   faces=[ [0,1,4],[1,2,4],[2,3,4],[3,0,4],              // each triangle side
#               [1,0,3],[2,1,3] ]                         // two triangles for square base
#  );
h = polygon(points=[(0,0,0),(5,0,0),(5,5,0),(0,5,0)])
h = h.linear_extrude(height=10)
h.save_as_scad()