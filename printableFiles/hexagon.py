import math
from solid2 import *
from typing import Sequence as _Sequence, \
    Tuple as _Tuple, \
    Union as _Union
import numpy as np
from solid2.extensions.bosl2 import cuboid, regular_ngon, \
    RIGHT, FRONT, BACK, BOT, CENTER, \
    path_sweep, xdistribute

from solid2.extensions.bosl2 import gears, beziers, screws, cubetruss

from solid2.extensions.bosl2 import gears, beziers, screws, cubetruss

##############
## Settings ##
##############

resolution = 100
hexagon_height = 2.0
hexagon_outer_width = 5.1961525

hexagon_bavel_size = 1.26

magnet_depth = 0.15  # inclusive of printer accuracy
magnet_radius = 0.53  # inclusive of printer accuracy
magnet_height_over_ground = 0.25

# holes for toppings
index_of_hexagons_for_toppings = [1, 2]

# Street
generate_street = True
street_entry_hexagon_indecies = [(0, 4)]
street_indent_height = 0.1
street_width_scalar = 1

# Water
generate_water = True
water_entry_hexagon_indecies = [(3, 5)]
water_indent_height = -0.58
water_width_scalar = 0.75

##########
## Code ##
##########

set_global_fn(resolution)
set_global_fa(resolution)
set_global_fs(resolution)


def generate_hexagon_vertecies(radius: float) -> list[_Tuple[float, float]]:
    vertices = []
    for i in range(6):
        angle_rad = math.radians(60 * i)  # 60 degrees for each vertex
        x = radius * math.cos(angle_rad)
        y = radius * math.sin(angle_rad)
        vertices.append((x, y))
    return vertices


def subdivide_hexagon(hexagon_vertices):
    center = (0, 0)  # Assuming the hexagon is centered at (0, 0)
    triangles = []

    # Create triangles by connecting center to each pair of adjacent vertices
    num_vertices = len(hexagon_vertices)
    for i in range(num_vertices):
        v1 = hexagon_vertices[i]
        # Wrap around to connect the last vertex to the first
        v2 = hexagon_vertices[(i + 1) % num_vertices]

        # Append triangle vertices (center, v1, v2) as flattened tuples
        triangles.append(center)
        triangles.append(v1)
        triangles.append(v2)

    return triangles


def generate_hexagon(radius: float, center: _Tuple[float, float] | None = None):
    # Calculate the coordinates of the vertices of a regular hexagon
    org_vertices = generate_hexagon_vertecies(radius=radius)
    subdivided_vertices = subdivide_hexagon(org_vertices)

    # Create a polygon using the calculated vertices
    hexagon = polygon(points=subdivided_vertices)
    if center is not None:
        hexagon = hexagon.translateX(center[0]).translateY(center[1])
        for i in range(org_vertices.__len__()):
            org_vertices[i] = (org_vertices[i][0] + center[0],
                               org_vertices[i][1] + center[1])
        for i in range(subdivided_vertices.__len__()):
            subdivided_vertices[i] = (subdivided_vertices[i][0] + center[0],
                                      subdivided_vertices[i][1] + center[1])

    return hexagon, org_vertices


def convertOuterHexagonSizeToInnerSize(outerHexagonSize: float):
    angle_rad = math.radians(60 * 0)
    temp_x_1 = outerHexagonSize * math.cos(angle_rad)
    temp_y_1 = outerHexagonSize * math.sin(angle_rad)
    angle_rad = math.radians(60 * 1)
    temp_x_2 = outerHexagonSize * math.cos(angle_rad)
    temp_y_2 = outerHexagonSize * math.sin(angle_rad)
    innerHexagonSizeX = temp_x_1 + 0.5 * (temp_x_2-temp_x_1)
    innerHexagonSizeY = temp_y_1 + 0.5 * (temp_y_2-temp_y_1)
    innerHexagonSize = math.sqrt(innerHexagonSizeX**2 + innerHexagonSizeY**2)
    return innerHexagonSize


innerHexagonSize = convertOuterHexagonSizeToInnerSize(hexagon_outer_width)

h = generate_hexagon(hexagon_outer_width)[0]
h = h.linear_extrude(height=hexagon_height)

all_hexagon_objects = [h]

outerHexagonVertecies: list[list[_Tuple[float, float]]] = []

for i in range(6):
    angle_rad = math.radians(60 * i + 30)  # 60 degrees for each vertex
    x = hexagon_outer_width * math.cos(angle_rad)
    y = hexagon_outer_width * math.sin(angle_rad)
    magnitude = 1.0/math.sqrt(x**2 + y**2)
    hexa, vert = generate_hexagon(hexagon_outer_width, center=(
        magnitude*x*innerHexagonSize*2, magnitude*y*innerHexagonSize*2))
    outerHexagonVertecies.append(vert)
    hexa = hexa.linear_extrude(height=hexagon_height)
    all_hexagon_objects.append(hexa)


# remove topping holes
subtractions_from_model = []
for hexagon_index_for_hole in index_of_hexagons_for_toppings:
    # construct high res cylinder for boolean operation
    removeTool = cylinder(h=magnet_depth*4,
                          center=True, r=magnet_radius)
    removeTool = removeTool.translateZ(hexagon_height)

    # find position of hole
    if (hexagon_index_for_hole == 0):
        position = (0, 0)
    else:
        # 60 degrees for each vertex
        angle_rad = math.radians(60 * hexagon_index_for_hole + 30)
        x = hexagon_outer_width * math.cos(angle_rad)
        y = hexagon_outer_width * math.sin(angle_rad)
        magnitude = 1.0/math.sqrt(x**2 + y**2)
        position = (magnitude*x*innerHexagonSize*2,
                    magnitude*y*innerHexagonSize*2)

    removeTool = removeTool.translateX(position[0]).translateY(position[1])
    subtractions_from_model.append(removeTool)

removeTool = union()(subtractions_from_model)

for i in range(all_hexagon_objects.__len__()):
    all_hexagon_objects[i] = all_hexagon_objects[i] - removeTool


def shift_line_3d(line, x, y, z):
    p1, p2 = line  # Unpack the line into two endpoints

    # Calculate direction vector from p1 to p2
    direction_vector = np.array([p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]])

    # Calculate magnitude of the direction vector
    magnitude = np.linalg.norm(direction_vector)

    if magnitude == 0:
        return line  # Handle zero-length line case

    # Normalize the direction vector to get a unit vector
    unit_direction = direction_vector / magnitude

    # Calculate two perpendicular vectors
    if unit_direction[0] == 0 and unit_direction[1] == 0:
        perpendicular_vector1 = np.array([1, 0, 0])
    else:
        perpendicular_vector1 = np.array(
            [-unit_direction[1], unit_direction[0], 0])
    perpendicular_vector1 /= np.linalg.norm(perpendicular_vector1)

    perpendicular_vector2 = np.cross(unit_direction, perpendicular_vector1)

    # Shift endpoints along the perpendicular vectors
    translated_p1 = p1 + x * perpendicular_vector1 + \
        y * perpendicular_vector2 + z * unit_direction
    translated_p2 = p2 + x * perpendicular_vector1 + \
        y * perpendicular_vector2 + z * unit_direction

    # Create the shifted line
    shifted_line = (tuple(translated_p1), tuple(translated_p2))

    return shifted_line


def addBevel(object, line: _Tuple[_Tuple[float, float, float], _Tuple[float, float, float]], depth: float):

    line = ((line[0][0] + 0.0001*(line[1][0]-line[0][0]), line[0][1] + 0.0001*(line[1][1]-line[0][1]), line[0][2] + 0.0001*(line[1][2]-line[0][2])),
            (line[1][0] + 0.0001*(line[0][0]-line[1][0]), line[1][1] + 0.0001*(line[0][1]-line[1][1]), line[1][2] + 0.0001*(line[0][2]-line[1][2])))

    shifted_line_1 = shift_line_3d(line, depth*0.75, 0, 0)
    shifted_line_2 = shift_line_3d(line, -depth*0.75, 0, 0)
    shifted_line_3 = shift_line_3d(line, 0, -depth, 0)

    points = [
        shifted_line_1[0],  # 0
        shifted_line_1[1],  # 1
        shifted_line_2[0],  # 2
        shifted_line_2[1],  # 3
        shifted_line_3[0],  # 4
        shifted_line_3[1],  # 5
    ]

    # first, construct a triangle on that line
    triangle = polyhedron(
        points=points,
        faces=[
            (0, 2, 4), (1, 5, 3),  # end caps
            (0, 4, 1), (4, 5, 1), (4, 2, 3), (5, 4, 3),  # side panels
            (2, 0, 1), (3, 2, 1)  # top face
        ])

    # then subtract triangle from object
    if object is None:
        return triangle.translateZ(hexagon_height)
    else:
        object = object-(triangle.translateZ(hexagon_height))
        return object


def angle_with_x_axis(lineForHole: _Tuple[_Tuple[float, float], _Tuple[float, float]]):
    point1, point2 = lineForHole
    x1, y1 = point1
    x2, y2 = point2

    # Calculate direction vector
    dx = x2 - x1
    dy = y2 - y1

    # Calculate angle with x-axis using atan2
    angle_rad = math.atan2(dy, dx)

    # Convert angle from radians to degrees (optional)
    angle_deg = math.degrees(angle_rad)
    return angle_rad, angle_deg


def addMagnetHoleOnSide(object, magnetHoleHeightAboveGround: float, magnetHoleRadius: float, magnetHoleDepth: float, lineForHole: _Tuple[_Tuple[float, float], _Tuple[float, float]]):
    # construct high res cylinder for boolean operation
    cylinderTool = cylinder(h=magnetHoleDepth*4,
                            center=True, r=magnetHoleRadius)

    centerOfLine = (lineForHole[1][0] - 0.5 * (lineForHole[1][0] - lineForHole[0][0]),
                    lineForHole[1][1] - 0.5 * (lineForHole[1][1] - lineForHole[0][1]))

    # rotate tool along line
    _, angle_deg = angle_with_x_axis(lineForHole)
    cylinderTool = cylinderTool.rotateX(90).rotateZ(angle_deg).translateX(centerOfLine[0]).translateY(
        centerOfLine[1]).translateZ(magnetHoleHeightAboveGround + magnetHoleRadius)

    return object - cylinderTool


tools_to_remove = []
tools_settings = []
for hexVerts in outerHexagonVertecies:
    length = hexVerts.__len__()
    for i in range(length):
        p1 = hexVerts[i]
        p2 = hexVerts[(i+1) % length]

        # small optimization to no add holes everywhere... still not perfect though
        angle_of_line = np.cross(
            np.array([p1[0], p1[1]]), np.array([p2[0], p2[1]]))
        if angle_of_line > 0.1:
            for i in range(all_hexagon_objects.__len__()):
                all_hexagon_objects[i] = addMagnetHoleOnSide(all_hexagon_objects[i], magnet_height_over_ground,
                                                             magnet_radius, magnet_depth, [p1, p2])

        # # TODO add me after everything is generated...
        bevel_tool_settings = ((p1[0], p1[1], 1.2), (p2[0], p2[1], 1.2))
        if bevel_tool_settings not in tools_settings:
            tool = addBevel(None, bevel_tool_settings, hexagon_bavel_size)
            tools_to_remove.append(tool)
            tools_settings.append(bevel_tool_settings)
        else:
            print("skipped")

for i in range(all_hexagon_objects.__len__()):
    all_hexagon_objects[i] = difference()(
        all_hexagon_objects[i], union()(tools_to_remove))


def getOuterHexFlowerLines(index: int) -> _Tuple[_Tuple[_Tuple[float, float], _Tuple[float, float]], _Tuple[_Tuple[float, float], _Tuple[float, float]], _Tuple[_Tuple[float, float], _Tuple[float, float]]]:
    # find matching hexagon indecies
    top_hexagon_index = index % 6
    hexagon_right_of_top_hex_index = (index - 1) % 6

    top_hexagon = outerHexagonVertecies[top_hexagon_index]
    hexagon_right_of_top_hex = outerHexagonVertecies[hexagon_right_of_top_hex_index]

    top_hex_line = (top_hexagon[(2 + index-1) % 6],
                    top_hexagon[(1 + index-1) % 6])
    top_right_hex_line = (
        top_hexagon[(2 + index-1-1) % 6], top_hexagon[(2 + index-1-2) % 6])
    adjecent_hexagon_right_of_top_hex_line = (hexagon_right_of_top_hex[(
        2 + index - 1) % 6], hexagon_right_of_top_hex[(1 + index - 1) % 6])

    return (top_hex_line, top_right_hex_line, adjecent_hexagon_right_of_top_hex_line)


def getCenterOfThreeLines(lineA: _Tuple[_Tuple[float, float], _Tuple[float, float]], lineB: _Tuple[_Tuple[float, float], _Tuple[float, float]], lineC: _Tuple[_Tuple[float, float], _Tuple[float, float]]) -> tuple[float, float]:
    center_of_three_lines_x = lineA[0][0] + 0.5 * (lineC[1][0]-lineA[0][0])
    center_of_three_lines_y = lineA[0][1] + 0.5 * (lineC[1][1]-lineA[0][1])
    center_of_three_lines = (center_of_three_lines_x, center_of_three_lines_y)
    return center_of_three_lines


if generate_street is True:
    subtractions_from_model = []
    counter = 0
    for streetIndexPair in street_entry_hexagon_indecies:
        counter += 1
        street_entry_hexagon_index, street_exit_hexagon_index = streetIndexPair
        entry_line_a, entry_line_b, entry_line_c = getOuterHexFlowerLines(
            street_entry_hexagon_index)
        center_of_entry_path = getCenterOfThreeLines(
            entry_line_a, entry_line_b, entry_line_c)

        exit_line_a, exit_line_b, exit_line_c = getOuterHexFlowerLines(
            street_exit_hexagon_index)
        center_of_exit_path = getCenterOfThreeLines(
            exit_line_a, exit_line_b, exit_line_c)

        entry_line_direction = np.array(
            [entry_line_c[1][0], entry_line_c[1][1], 0],) - np.array([entry_line_a[0][0], entry_line_a[0][1], 0])
        entry_line_direction_magnitude = np.linalg.norm(entry_line_direction)

        exit_line_direction = np.array(
            [exit_line_c[1][0], exit_line_c[1][1], 0],) - np.array([exit_line_a[0][0], exit_line_a[0][1], 0])
        exit_line_direction_magnitude = np.linalg.norm(exit_line_direction)

        unit_direction = np.array([0, 0, -1])
        perpendicular_vector_entry = np.cross(
            unit_direction, entry_line_direction * 1.0/entry_line_direction_magnitude)
        perpendicular_vector_exit = np.cross(
            unit_direction, exit_line_direction * 1.0/exit_line_direction_magnitude)

        inline_controll_point_a = np.array(
            [center_of_entry_path[0], center_of_entry_path[1], 0]) + 11 * perpendicular_vector_entry
        inline_controll_point_b = np.array(
            [center_of_exit_path[0], center_of_exit_path[1], 0]) + 11 * perpendicular_vector_exit

        exit_to_entry_vector = np.array([center_of_exit_path[0], center_of_exit_path[1], 0]) - np.array([
            center_of_entry_path[0], center_of_entry_path[1], 0])
        size_of_exit_to_entry_vector = np.linalg.norm(exit_to_entry_vector)

        center_of_hexagon = np.array([0., 0., 0.])

        midpoint_of_street = np.array([center_of_exit_path[0], center_of_exit_path[1], 0]) - \
            np.array([center_of_entry_path[0], center_of_entry_path[1], 0])
        midpoint_of_street *= 0.1 * 1/np.linalg.norm(midpoint_of_street)

        # Calculate magnitude of the direction vector
        street_width = np.linalg.norm(np.array([entry_line_c[1][0], entry_line_c[1][1]],) - np.array(
            [entry_line_a[0][0], entry_line_a[0][1]])) * street_width_scalar
        # since street shape is turned 45 degree, we must widen the street:
        diagonaled_street_width = math.sqrt(((street_width)**2)*2)

        # Add first street using curves???
        sbez = [
            [center_of_entry_path[0], center_of_entry_path[1]],
            [-midpoint_of_street[0], -midpoint_of_street[1]],
            [-center_of_hexagon[0], -center_of_hexagon[1]],
            [center_of_exit_path[0], center_of_exit_path[1]],
        ]

        street_tool = path_sweep(regular_ngon(n=4, d=diagonaled_street_width, spin=45), beziers.bezpath_curve(
            sbez, N=sbez.__len__()-1, splinesteps=resolution))

        street_tool = street_tool.recolor("#99f")

        street_tool_entry_points = [
            np.array(sbez[0])*0.9,
            np.array([0., 0.]) + np.array([center_of_entry_path[0],
                                           center_of_entry_path[1]]) * 5,
            np.array([0., 0.]) + np.array([center_of_entry_path[0],
                                           center_of_entry_path[1]]) * 10,
        ]
        street_tool_entry = path_sweep(regular_ngon(n=4, d=diagonaled_street_width, spin=45), beziers.bezpath_curve(
            street_tool_entry_points, N=street_tool_entry_points.__len__()-1, splinesteps=resolution))
        street_tool += street_tool_entry

        street_tool_ending_points = [
            np.array(sbez[sbez.__len__()-1])*0.9,
            np.array([0., 0.]) + np.array([center_of_exit_path[0],
                                           center_of_exit_path[1]]) * 5,
            np.array([0., 0.]) + np.array([center_of_exit_path[0],
                                           center_of_exit_path[1]]) * 10,
        ]
        street_tool_ending = path_sweep(regular_ngon(n=4, d=diagonaled_street_width, spin=45), beziers.bezpath_curve(
            street_tool_ending_points, N=street_tool_ending_points.__len__()-1, splinesteps=resolution))
        street_tool += street_tool_ending

        # add collection of bevels
        bevelsTool = None
        for hexVerts in outerHexagonVertecies:
            length = hexVerts.__len__()
            for i in range(length):
                p1 = hexVerts[i]
                p2 = hexVerts[(i+1) % length]

                if bevelsTool is None:
                    bevelsTool = addBevel(
                        None, ((p1[0], p1[1], 1.2), (p2[0], p2[1], 1.2)), hexagon_bavel_size)
                else:
                    bevelsTool += addBevel(
                        None, ((p1[0], p1[1], 1.2), (p2[0], p2[1], 1.2)), hexagon_bavel_size)

        # shrink bevels to street box sizes(union???)
        bevelsTool = intersection()(bevelsTool, street_tool)

        # combine bevels with street box (with correct translation)
        street_tool = street_tool.translateZ(
            street_width/2 + hexagon_height - street_indent_height + counter*0.00001)
        subtractions_from_model.append(street_tool)

        bevelsTool = bevelsTool.translateZ(-street_indent_height)
        subtractions_from_model.append(bevelsTool)

    removeTool = union()(subtractions_from_model)

    for i in range(all_hexagon_objects.__len__()):
        all_hexagon_objects[i] = all_hexagon_objects[i] - removeTool


if generate_water is True:
    subtractions_from_model = []
    counter = 0
    for waterIndexPair in water_entry_hexagon_indecies:
        counter += 1
        water_entry_hexagon_index, water_exit_hexagon_index = waterIndexPair
        entry_line_a, entry_line_b, entry_line_c = getOuterHexFlowerLines(
            water_entry_hexagon_index)
        center_of_entry_path = getCenterOfThreeLines(
            entry_line_a, entry_line_b, entry_line_c)

        exit_line_a, exit_line_b, exit_line_c = getOuterHexFlowerLines(
            water_exit_hexagon_index)
        center_of_exit_path = getCenterOfThreeLines(
            exit_line_a, exit_line_b, exit_line_c)

        entry_line_direction = np.array(
            [entry_line_c[1][0], entry_line_c[1][1], 0],) - np.array([entry_line_a[0][0], entry_line_a[0][1], 0])
        entry_line_direction_magnitude = np.linalg.norm(entry_line_direction)

        exit_line_direction = np.array(
            [exit_line_c[1][0], exit_line_c[1][1], 0],) - np.array([exit_line_a[0][0], exit_line_a[0][1], 0])
        exit_line_direction_magnitude = np.linalg.norm(exit_line_direction)

        unit_direction = np.array([0, 0, -1])
        perpendicular_vector_entry = np.cross(
            unit_direction, entry_line_direction * 1.0/entry_line_direction_magnitude)
        perpendicular_vector_exit = np.cross(
            unit_direction, exit_line_direction * 1.0/exit_line_direction_magnitude)

        inline_controll_point_a = np.array(
            [center_of_entry_path[0], center_of_entry_path[1], 0]) + 11 * perpendicular_vector_entry
        inline_controll_point_b = np.array(
            [center_of_exit_path[0], center_of_exit_path[1], 0]) + 11 * perpendicular_vector_exit

        exit_to_entry_vector = np.array([center_of_exit_path[0], center_of_exit_path[1], 0]) - np.array([
            center_of_entry_path[0], center_of_entry_path[1], 0])
        size_of_exit_to_entry_vector = np.linalg.norm(exit_to_entry_vector)

        center_of_hexagon = np.array([0., 0., 0.])

        midpoint_of_water = np.array([center_of_exit_path[0], center_of_exit_path[1], 0]) - \
            np.array([center_of_entry_path[0], center_of_entry_path[1], 0])
        midpoint_of_water *= 0.1 * 1/np.linalg.norm(midpoint_of_water)

        # Calculate magnitude of the direction vector
        water_width = np.linalg.norm(np.array([entry_line_c[1][0], entry_line_c[1][1]],) - np.array(
            [entry_line_a[0][0], entry_line_a[0][1]])) * water_width_scalar
        # since water shape is turned 45 degree, we must widen the water:
        diagonaled_water_width = math.sqrt(((water_width)**2)*2)

        # Add first water using curves???
        sbez = [
            [center_of_entry_path[0], center_of_entry_path[1]],
            [-midpoint_of_water[0], -midpoint_of_water[1]],
            [-center_of_hexagon[0], -center_of_hexagon[1]],
            [center_of_exit_path[0], center_of_exit_path[1]],
        ]

        water_tool = path_sweep(regular_ngon(n=6, d=diagonaled_water_width, spin=60), beziers.bezpath_curve(
            sbez, N=sbez.__len__()-1, splinesteps=resolution))

        water_tool = water_tool.recolor("#99f")

        water_tool_entry_points = [
            np.array(sbez[0])*0.9,
            np.array([0., 0.]) + np.array([center_of_entry_path[0],
                                           center_of_entry_path[1]]) * 5,
            np.array([0., 0.]) + np.array([center_of_entry_path[0],
                                           center_of_entry_path[1]]) * 10,
        ]
        water_tool_entry = path_sweep(regular_ngon(n=6, d=diagonaled_water_width, spin=60), beziers.bezpath_curve(
            water_tool_entry_points, N=water_tool_entry_points.__len__()-1, splinesteps=resolution))
        water_tool += water_tool_entry

        water_tool_ending_points = [
            np.array(sbez[sbez.__len__()-1])*0.9,
            np.array([0., 0.]) + np.array([center_of_exit_path[0],
                                           center_of_exit_path[1]]) * 5,
            np.array([0., 0.]) + np.array([center_of_exit_path[0],
                                           center_of_exit_path[1]]) * 10,
        ]
        water_tool_ending = path_sweep(regular_ngon(n=6, d=diagonaled_water_width, spin=60), beziers.bezpath_curve(
            water_tool_ending_points, N=water_tool_ending_points.__len__()-1, splinesteps=resolution))
        water_tool += water_tool_ending

        # add collection of bevels
        bevelsTool = None
        for hexVerts in outerHexagonVertecies:
            length = hexVerts.__len__()
            for i in range(length):
                p1 = hexVerts[i]
                p2 = hexVerts[(i+1) % length]

                if bevelsTool is None:
                    bevelsTool = addBevel(
                        None, ((p1[0], p1[1], 1.2), (p2[0], p2[1], 1.2)), hexagon_bavel_size)
                else:
                    bevelsTool += addBevel(
                        None, ((p1[0], p1[1], 1.2), (p2[0], p2[1], 1.2)), hexagon_bavel_size)

        # shrink bevels to water box sizes(union???)
        bevelsTool = intersection()(bevelsTool, water_tool)

        # combine bevels with water box (with correct translation)
        water_tool = water_tool.translateZ(
            water_width/2 + hexagon_height - water_indent_height + counter*0.00001)
        subtractions_from_model.append(water_tool)

        bevelsTool = bevelsTool.translateZ(-water_indent_height)
        # subtractions_from_model.append(bevelsTool)

    removeTool = union()(subtractions_from_model)

    for i in range(all_hexagon_objects.__len__()):
        all_hexagon_objects[i] = all_hexagon_objects[i] - removeTool


h = union()(all_hexagon_objects)

print("finished descibing model")
h.scale(5).save_as_scad()
print("finished save_as_scad")
# takes ages and does not render bezier stuff...
# h = h.scale(5)
# h.save_as_stl()
print("Ran")
