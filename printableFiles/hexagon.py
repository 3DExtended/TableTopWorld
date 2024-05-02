import math
from solid2 import *
from typing import Sequence as _Sequence,\
                   Tuple as _Tuple,\
                   Union as _Union
import numpy as np

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
hex_org = generate_hexagon(hexagonSize)

for i in range(6):
    angle_rad = math.radians(60 * i + 30)  # 60 degrees for each vertex
    x = hexagonSize * math.cos(angle_rad)
    y = hexagonSize * math.sin(angle_rad)
    magnitude = 1.0/math.sqrt(x**2 + y**2)
    
    h = h + generate_hexagon(hexagonSize, center=(magnitude*x*innerHexagonSize*2,magnitude*y*innerHexagonSize*2))


def shift_line_3d(line, x, y, z):
    p1, p2 = line  # Unpack the line into two endpoints
    
    # Calculate direction vector from p1 to p2
    direction_vector = np.array([p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]])
    
    # Calculate magnitude of the direction vector
    magnitude = np.linalg.norm(direction_vector)
    
    if magnitude == 0:
        print("magnitude was 0")
        return line  # Handle zero-length line case
    
    # Normalize the direction vector to get a unit vector
    unit_direction = direction_vector / magnitude
    
    # Calculate two perpendicular vectors
    if unit_direction[0] == 0 and unit_direction[1] == 0:
        perpendicular_vector1 = np.array([1, 0, 0])
        print("a")
    else:
        perpendicular_vector1 = np.array([-unit_direction[1], unit_direction[0], 0])
    perpendicular_vector1 /= np.linalg.norm(perpendicular_vector1)
    
    perpendicular_vector2 = np.cross(unit_direction, perpendicular_vector1)
    print("perpendicular_vector2")
    print(perpendicular_vector2)
    
    # Shift endpoints along the perpendicular vectors
    translated_p1 = p1 + x * perpendicular_vector1 + y * perpendicular_vector2 + z * unit_direction
    translated_p2 = p2 + x * perpendicular_vector1 + y * perpendicular_vector2 + z * unit_direction
    
    # Create the shifted line
    shifted_line = (tuple(translated_p1), tuple(translated_p2))
    
    return shifted_line



def addBevel(object, line: _Tuple[_Tuple[float, float, float], _Tuple[float, float, float]], depth: float):
    # Example usage:
    shifted_line_1 = shift_line_3d(line, depth*0.75, 0, 0)
    shifted_line_2 = shift_line_3d(line, -depth*0.75, 0, 0)
    shifted_line_3 = shift_line_3d(line, 0, -depth, 0)
    
    print(line)
    print(shifted_line_3)
    
    # first, construct a triangle on that line
    triangle = polyhedron(
        points=[
            shifted_line_1[0], # 0
            shifted_line_1[1], # 1
            shifted_line_2[0], # 2
            shifted_line_2[1], # 3
            shifted_line_3[0], # 4
            shifted_line_3[1], # 5
        ],
        faces=[
            (0,2,4), (1,5,3), # end caps
            (0,4,1),(4,5,1),(4,2,3),(5,4,3), # side panels
            (2,0,1), (3,2,1) # top face
        ])
    
    # then subtract triangle from object
    object = object-triangle
    return object

h = h.linear_extrude(height=1)

h = addBevel(h, ((0,0,1.2),(5,0,1.2)), 0.3)

h.save_as_scad()
h.save_as_stl()

print("Ran")