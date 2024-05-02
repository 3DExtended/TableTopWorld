import math
from solid2 import *
from typing import Sequence as _Sequence,\
                   Tuple as _Tuple,\
                   Union as _Union
                   

set_global_fn(100)

# d = cube(5) + sphere(5).right(5) - cylinder(r=2, h=6)
# polyhedron(
#   points=[ [10,10,0],[10,-10,0],[-10,-10,0],[-10,10,0], // the four points at base
#            [0,0,10]  ],                                 // the apex point 
#   faces=[ [0,1,4],[1,2,4],[2,3,4],[3,0,4],              // each triangle side
#               [1,0,3],[2,1,3] ]                         // two triangles for square base
#  );
# h = polygon(points=[(0,0,0),(5,0,0),(5,5,0),(0,5,0)])

def generate_hexagon_vertecies(radius:float) -> list[_Tuple[float, float]]:
    vertices = []
    for i in range(6):
        angle_rad = math.radians(60 * i)  # 60 degrees for each vertex
        x = radius * math.cos(angle_rad)
        y = radius * math.sin(angle_rad)
        vertices.append((x, y))
    return vertices

def generate_hexagon(radius:float, center:_Tuple[float, float]|None=None):
    # Calculate the coordinates of the vertices of a regular hexagon
    vertices = generate_hexagon_vertecies(radius=radius)
    
    # Create a polygon using the calculated vertices
    hexagon = polygon(points=vertices)
    if center is not None:
        hexagon = hexagon.translateX(center[0]).translateY(center[1])
    return hexagon

# furthest points of hexagon
hexagonSize = 5.1961525

def convertOuterHexagonSizeToInnerSize(outerHexagonSize: float):
    angle_rad = math.radians(60 * 0)
    temp_x_1 = hexagonSize * math.cos(angle_rad)
    temp_y_1 = hexagonSize * math.sin(angle_rad)
    angle_rad = math.radians(60 * 1)
    temp_x_2 = hexagonSize * math.cos(angle_rad)
    temp_y_2 = hexagonSize * math.sin(angle_rad)
    innerHexagonSizeX = temp_x_1 + 0.5 * (temp_x_2-temp_x_1)
    innerHexagonSizeY = temp_y_1 + 0.5 * (temp_y_2-temp_y_1)
    innerHexagonSize = math.sqrt(innerHexagonSizeX**2 + innerHexagonSizeY**2)
    return innerHexagonSize

innerHexagonSize = convertOuterHexagonSizeToInnerSize(hexagonSize)

h = generate_hexagon(hexagonSize)

for i in range(6):
    angle_rad = math.radians(60 * i + 30)  # 60 degrees for each vertex
    x = hexagonSize * math.cos(angle_rad)
    y = hexagonSize * math.sin(angle_rad)
    magnitude = 1.0/math.sqrt(x**2 + y**2)
    
    h = h + generate_hexagon(hexagonSize, center=(magnitude*x*innerHexagonSize*2,magnitude*y*innerHexagonSize*2))

h = h.linear_extrude(height=1)

h.save_as_scad()