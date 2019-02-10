# coding: utf-8
import os 
import math
import pymol

import numpy as np

from Bio.PDB.PDBParser import PDBParser

import bio_project.utils as utils


def center_mass_unit(draw, unit_dir):
    """
    compute center of mass for each unit

    :param draw: boolean, if True pymol is used to draw the respective 
                 geometric element
    :return: dictionary containing coordinates of the centers of mass
    """
    unit_dir_list = os.listdir(unit_dir)
    unit_dir_list.sort(key=utils.natural_keys) 
    centers = {}  
    for unit in unit_dir_list:
        pdb_parser = PDBParser()
        structure = pdb_parser.get_structure(unit, unit_dir + unit)
        # computing center of mass 
        tot_mass = 0
        x = 0
        y = 0
        z = 0
        for atom in structure.get_atoms():
            x += atom.mass*atom.coord[0]
            y += atom.mass*atom.coord[1]
            z += atom.mass*atom.coord[2]
            tot_mass += atom.mass

        x = x/tot_mass
        y = y/tot_mass
        z = z/tot_mass

        centers[unit] = (x, y, z)

        # drawing center of mass
        if draw:
            utils.draw_center_of_mass(centers)

    return centers


def distance_center_of_mass(centers, draw, unit_dir):
    """
    compute distance between each center of mass

    :param centers: list of coordinates of the centers of mass for each unit
    :param draw: boolean, if True pymol is used to draw the respective 
                 geometric element
    :return: dictionary containig distances between center of mass
    """
    unit_dir_list = os.listdir(unit_dir)
    unit_dir_list.sort(key=utils.natural_keys) 
    distances = {}
    for i in range(len(centers)):
        unit_1 = unit_dir_list[i]
        distances[unit_1] = {}
        for j in range(len(centers)):
            if i != j:
                unit_2 = unit_dir_list[j]
                distances[unit_1][unit_2] = distance(centers[unit_1],
                                                     centers[unit_2])

    # drawing distances
    if draw:
        utils.draw_distance_center_of_mass(centers)

    return distances


def distance_alpha_c(centers, draw, unit_dir):
    """
    compute the distance between center of mass and each alpha carbon 
    in the corresponding unit 

    :param centers: list of coordinates of the centers of mass for each unit
    :param draw: boolean, if True pymol is used to draw the respective 
                 geometric element
    :return: dictionary of the distaces between center of mass and alpha carbon
    """
    unit_dir_list = os.listdir(unit_dir)
    unit_dir_list.sort(key=utils.natural_keys) 
    distances = {}
    for unit in unit_dir_list:
        pdb_parser = PDBParser()
        structure = pdb_parser.get_structure(unit, unit_dir + unit)
        alpha_c = alpha_carbon(structure)
        center = centers[unit]
        distances[unit] = {}
        for i, ca in enumerate(alpha_c):
            distances[unit][ca.get_id() + str(i)] = distance(ca.get_coord(),
                                                             center)

        # drawing distances
        if draw:
            utils.draw_distance_center_mass_alpha(unit, center, alpha_c)

    return distances


def handedness(centers, draw, unit_dir):
    """
    compute the handeness (Left or Right) of each unit 

    :param centers: list of coordinates of the centers of mass for each unit
    :param draw: boolean, if True pymol is used to draw the respective 
                 geometric element
    :return: dictionary containing the handedness for each unit
    """
    unit_dir_list = os.listdir(unit_dir)
    unit_dir_list.sort(key=utils.natural_keys) 
    handedness = {}
    for i in range(len(centers)):
        unit = unit_dir_list[i]
        pdb_parser = PDBParser()
        structure = pdb_parser.get_structure(unit, unit_dir + unit)
        alpha_c = alpha_carbon(structure)
        ca_1 = alpha_c[0]
        ca_2 = alpha_c[1]
        x = np.array(ca_1.get_coord()) - np.array(centers[unit])
        v = np.array(ca_2.get_coord()) - np.array(centers[unit])
        y = np.array(centers[unit])
        if i != len(centers)-1:
            y = y + np.array(centers[unit_dir_list[i+1]]) 
        if i == len(centers)-1:
            y = 1 - (y + np.array(centers[unit_dir_list[i-1]]))
        z = np.cross(y, x)
        x_trans = np.linalg.solve(np.array([x, y, z]).T, x)
        v_trans = np.linalg.solve(np.array([x, y, z]).T, v)
        hand = np.cross(x_trans, v_trans)[2]
        
        if hand >= 0:
            handedness[unit] = "R"
        else:
            handedness[unit] = "L"

        # drawing vector used to compute handedness
        if draw:
            utils.draw_vector(list(ca_1.get_coord()), list(centers[unit]))
            utils.draw_vector(list(ca_2.get_coord()), list(centers[unit]))

    return handedness


def twist(centers, draw, unit_dir):
    """
    compute the rotation between units

    :param centers: list of coordinates of the centers of mass for each unit
    :param draw: boolean, if True pymol is used to draw the respective 
                 geometric element
    :return: dictionary containing twist between units
    """
    unit_dir_list = os.listdir(unit_dir)
    unit_dir_list.sort(key=utils.natural_keys) 
    twist = {}
    for i in range(len(centers)-1):
        unit_1 = unit_dir_list[i]
        unit_2 = unit_dir_list[i+1]
        pdb_parser = PDBParser()
        # get coordinates first unit
        structure_1 = pdb_parser.get_structure(unit_1, unit_dir + unit_1)
        center_1 = np.array(centers[unit_1])
        ca_1 = np.array(alpha_carbon(structure_1)[0].get_coord())
        # get coordinates second unit
        structure_2 = pdb_parser.get_structure(unit_2, unit_dir + unit_2)
        center_2 = np.array(centers[unit_2])
        ca_2 = np.array(alpha_carbon(structure_2)[0].get_coord())
        # build vector representing the plane approximating the geometry
        # of the unit
        v_1 = ca_1 - center_1
        v_2 = ca_2 - center_2
        # compute the angle between the two planes
        prod_sum_wise = np.absolute(np.sum(v_1 * v_2))
        sqrt_1 = np.sqrt(np.sum(v_1**2))
        sqrt_2 = np.sqrt(np.sum(v_2**2))
        angle = prod_sum_wise/(sqrt_1*sqrt_2)
        np.degrees(np.arccos(angle))
        twist[unit_1 + '_' + unit_2] = angle

        #  drawing vector used to compute twist
        if draw:
            utils.draw_vector(list(ca_1), list(ca_2))
            utils.draw_vector(list(center_1), list(ca_1))
            utils.draw_vector(list(center_2), list(ca_2))

    return twist


def pitch(centers, draw, unit_dir):
    """
    compute the pitch of the protein, vertical angle between vectors
    connecting units

    :param centers: list of coordinates of the centers of mass for each unit
    :param draw: boolean, if True pymol is used to draw the respective
                 geometric element
    :return: dictionary containing pitch
    """
    unit_dir_list = os.listdir(unit_dir)
    unit_dir_list.sort(key=utils.natural_keys) 
    pitch = {}
    for i in range(len(centers)-2):
        v_1, v_2, v_3 = aux(i, centers, draw, unit_dir)
        pitch[unit_dir_list[i] + '_' + unit_dir_list[i+1]] = (angle(v_1, v_3) - 
                                                              angle(v_1, v_2))

    return pitch


def curvature(centers, draw, unit_dir):
    """
    compute the curvature of the protein, horizontal angle between vectors
    connecting units

    :param centers: list of coordinates of the centers of mass for each unit
    :param draw: boolean, if True pymol is used to draw the respective
                 geometric element
    :return: dictionary containing curvature
    """
    unit_dir_list = os.listdir(unit_dir)
    unit_dir_list.sort(key=utils.natural_keys) 
    curvature = {}
    for i in range(len(centers)-2):
        v_1, v_2, v_3 = aux(i, centers, draw, unit_dir)
        curvature[unit_dir_list[i] +
                  '_' + unit_dir_list[i+1]] = angle(v_2, v_3)

    return curvature


# -----------------------------------------------------------------------------
# AUXILIARY FUNCTIONS 
# -----------------------------------------------------------------------------


def distance(coord_1, coord_2):
    """
    compute the euclidean distance between two points

    :param coord_1: coordinates first point 
    :param coord_2: coordinates second point
    :return: distance between two points
    """
    return math.sqrt((coord_1[0] - coord_2[0])**2 +
                     (coord_1[1] - coord_2[1])**2 +
                     (coord_1[2] - coord_2[2])**2)


def angle(x, y):
    """
    compute the angle (in degrees) between two vectors 

    :param x: first vector
    :param y: second vector
    :return: angle in degrees between two vectors
    """
    dot = np.dot(x, y)
    norm_x = np.linalg.norm(x)
    norm_y = np.linalg.norm(y)
    temp_angle = dot/(norm_x*norm_y)
    return np.degrees(np.arccos(temp_angle))


def alpha_carbon(structure):
    """
    retrieve all the alpha carbon from the unit structure

    :param structure: bio python unit structure
    :return: list of alpha carbon cointained in the structure
    """
    return [atom for atom in structure.get_atoms() if atom.get_name() == 'CA']


def aux(i, centers, draw, unit_dir):
    """
    axiliary function, transforms vectors used to compute pitch and curvature,
    in order to have comparable arrays

    :param i: number used to identify the set of unit used
    :param centers: list of coordinates of the centers of mass for each unit
    :param draw: boolean, if True pymol is used to draw the respective
                 geometric element
    :return: 3 vector transposed according to a coordinate system that allows
             them to be compared
    """
    unit_dir_list = os.listdir(unit_dir)
    unit_dir_list.sort(key=utils.natural_keys) 
    unit_1 = unit_dir_list[i]
    unit_2 = unit_dir_list[i+1]
    unit_3 = unit_dir_list[i+2]
    pdb_parser = PDBParser()
    structure_1 = pdb_parser.get_structure(unit_1, unit_dir + unit_1)
    center_1 = np.array(centers[unit_1])
    ca_1 = np.array(alpha_carbon(structure_1)[-1].get_coord())  
    structure_2 = pdb_parser.get_structure(unit_2, unit_dir + unit_2)
    center_2 = np.array(centers[unit_2])
    structure_3 = pdb_parser.get_structure(unit_3, unit_dir + unit_3)
    center_3 = np.array(centers[unit_3])

    x = ca_1 - center_1
    y = center_2 - center_1
    v = center_3 - center_2
    z = np.cross(x, y)

    x_trans = np.linalg.solve(np.array([x, y, z]), x)
    y_trans = np.linalg.solve(np.array([x, y, z]), y)
    v_trans = np.linalg.solve(np.array([x, y, z]), v)

    if draw:
        utils.draw_vector(list(center_1), list(ca_1))
        utils.draw_vector(list(center_1), list(center_2))
        utils.draw_vector(list(center_2), list(center_3))

    return x_trans, y_trans, v_trans