# coding: utf-8
import re
import os 
import sys
import math

import pymol

import numpy as np

from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
from bio_project.utils import UnitSelection, parse_rdb, draw_center_of_mass, \
    draw_distance_center_of_mass, draw_distance_center_mass_alpha, draw_vector

pdb_f = "data/2z7x.pdb"
rdb_f = "data/rdb/2z7xb.db"
unit_dir = "data/unit/"

# parsing RDB protein file 
parsed_rdb_protein = parse_rdb(rdb_f)
# selecting and saving units of the rdb protein 
rdb_object = parse_rdb(rdb_f)
p = PDBParser()
pdb_structure = p.get_structure("PDB_structure", pdb_f)
io = PDBIO()
io.set_structure(pdb_structure[0][rdb_object['chain']])
for unit in parsed_rdb_protein['units']:
    range_unit = str(unit[0][1]) + '-' + str(unit[1][1])
    io.save(unit_dir + "unit" + '_' + range_unit + ".pdb",
            UnitSelection(unit[0][1], unit[1][1]))
# pymol
pymol.finish_launching()
for unit in os.listdir(unit_dir):
    pymol.cmd.load(unit_dir + unit, unit)


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [atoi(c) for c in re.split('(\d+)',text)]


def center_mass_unit(draw):
    """
    Compute center of mass for each unit
    """
    centers = {}
    for unit in os.listdir(unit_dir):
        pdb_parser = PDBParser()
        structure = pdb_parser.get_structure(unit, unit_dir + unit)
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
    
        if draw:
            draw_center_of_mass(centers)

    return centers


def distance(coord_1, coord_2):
    return math.sqrt((coord_1[0] - coord_2[0])**2 +
                     (coord_1[1] - coord_2[1])**2 +
                     (coord_1[2] - coord_2[2])**2)


def distance_center_of_mass(centers, draw):
    distances = []
    for unit_1 in centers:
        dist = []
        for unit_2 in centers:
            if unit_1 != unit_2:
                d = distance(centers[unit_1], centers[unit_2])
                dist.append(d)
        distances.append(dist)

    if draw:
        draw_distance_center_of_mass(centers)
    return distances


def distance_alpha_c(centers, draw):
    distances = []
    for unit in centers:
        pdb_parser = PDBParser()
        structure = pdb_parser.get_structure(unit, 
                                             unit_dir + unit)
        alpha_c = []
        for atom in structure.get_atoms():
            if atom.get_name() == 'CA':
                alpha_c.append(atom)
        dist = []
        for c in alpha_c:
            c_coord = c.get_coord()
            center = centers[unit]
            d = distance(c.get_coord(), center)
            dist.append(d) 
        distances.append(dist)

        if draw:
            draw_distance_center_mass_alpha(unit, center, alpha_c)

    return distances


def handedness(centers, draw):
    file_list = os.listdir(unit_dir)
    handedness_dic = {}
    for i, unit in enumerate(centers):
        pdb_parser = PDBParser()
        structure = pdb_parser.get_structure(unit, unit_dir + unit)
        ca_1 = None
        ca_2 = None
        i = 0 
        for atom in structure.get_atoms():
            if atom.get_name() == 'CA':
                if i == 0:
                    ca_1 = atom
                    i += 1
                else:
                    ca_2 = atom
                    break
        v_1 = np.array(ca_1.get_coord()) - np.array(centers[unit])
        draw_vector(list(ca_1.get_coord()), list(centers[unit]))
        v_2 = np.array(ca_2.get_coord()) - np.array(centers[unit])
        draw_vector(list(ca_2.get_coord()), list(centers[unit]))
        # z = None
        # if i == len(centers)-1:
        #     head = np.array(centers[unit])
        #     tail = np.array(centers[unit])
        #     z = 1 - (tail - head)
        # else:
        #     head = np.array(centers[file_list[i]])
        #     tail = np.array(centers[file_list[i+1]])
        #     z = tail - head
        # y = np.cross(z, v_1)

        # coord_system = np.array([v_1, y, z]).T
        # v_1_trans = np.linalg.solve(coord_system, v_1)
        # v_2_trans = np.linalg.solve(coord_system, v_2)
        cr_product_trans = np.cross(v_1, v_2)
        handedness = ""
        if cr_product_trans[2] >= 0:
            handedness = "R"
        else:
            handedness = "L"
        handedness_dic[unit] = handedness
    return handedness_dic


def twist(centers, draw):
    rotations = []
    list_dir = os.listdir(unit_dir)
    list_dir.sort(key=natural_keys)
    for i in range(len(centers)-1):
        unit_1 = list_dir[i]
        unit_2 = list_dir[i+1]
        pdb_parser = PDBParser()
        structure_1 = pdb_parser.get_structure(unit_1, unit_dir + unit_1)
        center_1 = np.array(centers[unit_1])
        ca_1 = None 
        for atom in structure_1.get_atoms():
            if atom.get_name() == "CA":
                ca_1 = np.array(atom.get_coord())
                break 
        draw_vector(list(center_1), list(ca_1))
        structure_2 = pdb_parser.get_structure(unit_2, unit_dir + unit_2)
        center_2 = np.array(centers[unit_2])
        ca_2 = None
        for atom in structure_2.get_atoms():
            if atom.get_name() == "CA":
                ca_2 = np.array(atom.get_coord())
                break 
        # c_center = np.cross(center_1, center_2)
        # c_ca = np.cross(ca_1, ca_2)
        # np.linalg.norm(c_center)
        # np.linalg.norm(c_ca)
        # rot_x = math.atan2(c_center[1], c_center[2])
        # rot_y = math.atan2(
        #     math.sqrt(c_ca[0]**2 + c_ca[1]**2), c_ca[2]) - \
        #     math.atan2(c_center[0], c_center[1])
        # rot_z = math.atan2(c_ca[1], c_ca[0]) - rot_x
        # rotations.append((rot_x, rot_y, rot_z))
        v_1 = ca_1 - center_1
        v_2 = ca_2 - center_2
        prod_sum_wise = np.sum(np.prod((v_1, v_2), axis=1))
        np.absolute(prod_sum_wise)
        sqrt_1 = np.sqrt(np.sum(v_1**2))
        sqrt_2 = np.sqrt(np.sum(v_2**2))
        temp_angle = prod_sum_wise/(sqrt_1 + sqrt_2)

        if draw:
            draw_vector(list(ca_1), list(ca_2))
            draw_vector(list(center_1), list(ca_1))
            draw_vector(list(center_2), list(ca_2))

    return np.degrees(np.arccos(temp_angle))
                

def angle(x, y):
    dot = np.dot(x, y)
    mag_x = np.linalg.norm(x)
    mag_y = np.linalg.norm(y)
    temp_angle = dot/(mag_x*mag_y)
    return np.degrees(np.arccos(temp_angle))


def curvature_pitch(centers, draw):
    curvature = []
    pitch = [] 
    list_dir = os.listdir(unit_dir)
    list_dir.sort(key=natural_keys)
    for i in range(len(centers)-2):
        unit_1 = list_dir[i]
        unit_2 = list_dir[i+1]
        unit_3 = list_dir[i+2]
        pdb_parser = PDBParser()
        structure_1 = pdb_parser.get_structure(unit_1, unit_dir + unit_1)
        center_1 = np.array(centers[unit_1])
        ca_1 = None 
        for atom in structure_1.get_atoms():
            if atom.get_name() == "CA":
                ca_1 = np.array(atom.get_coord())
                break 
        draw_vector(list(center_1), list(ca_1))
        structure_2 = pdb_parser.get_structure(unit_2, unit_dir + unit_2)
        center_2 = np.array(centers[unit_2])
        structure_3 = pdb_parser.get_structure(unit_3, unit_dir + unit_3)
        center_3 = np.array(centers[unit_3])

        v_1 = ca_1 - center_1
        v_2 = center_2 - center_1
        v_3 = center_3 - center_2

        curvature_ = angle(v_2, v_3)
        pitch_ = angle(v_1, v_3) - angle(v_1, v_2)

        curvature.append(curvature_)
        pitch.append(pitch_)
    
        if draw:
            draw_vector(list(center_1), list(ca_1))
            draw_vector(list(center_1), list(center_2))
            draw_vector(list(center_2), list(center_3))

    return curvature, pitch