# coding: utf-8
import re
# import pymol

from uuid import uuid4

from Bio.PDB import Select

# from pymol import cmd
# from pymol.cgo import *


class UnitSelection(Select):
    """
    class to handle the division of protein in unit
    """
    def __init__(self, start_id, end_id):
        """
        Initialization of the class

        :param start_id: identifier of the first element in the unit
        :param end_id: identifier of the last element in the unit
        """
        self.start_id = start_id
        self.end_id = end_id

    def accept_residue(self, residue):
        """
        Accept residue if in the range of initialization 

        :param residure: elment to check 
        :return: True if residue is in the interval [start_id, end_id]
        """
        hetflag, resseq, icode = residue.get_id()
        if self.start_id <= resseq <= self.end_id:
            return True
        return False


def atoi(text):
    """
    transform the a digit text into integers

    :param text: a string of text
    :return: an integer if the text is composed by digits,
             the text itself otherwise
    """
    return int(text) if text.isdigit() else text


def natural_keys(text):
    """
    split the text between character and digits

    :param text: a string of text
    :return: a list containig the splitted text
    """
    return [atoi(c) for c in re.split('(\d+)', text) if c != '']


def parse_rdb(filename):
    """
    parse RepeadDB protein file 

    :param filename: path to pdb file
    :return: dictionary representing the parsed file 
    """
    obj = {}
    with open(filename) as f:
        for line in f:
            line = line.strip().split()
            if line[0] == "SOURCE":
                obj['source'] = line[1]
            elif line[0] == "PDB":
                obj['id'] = line[1]
            elif line[0] == "CHAIN":
                obj['chain'] = line[1]
            elif line[0] == "REG":
                obj.setdefault('regions', [])
                obj['regions'].append((line[1], line[2], 
                                       line[3] + '.' + line[4]))
            elif line[0] == "UNIT":
                start, end = line[1:3]
                try:
                    start_id = (' ', int(start), ' ')
                except:
                    start_id = (' ', int(start[:-1]), start[-1])
                try:
                    end_id = (' ', int(end), ' ')
                except:
                    end_id = (' ', int(end[:-1]), end[-1])

                obj.setdefault('units', [])
                obj['units'].append((start_id, end_id))

            elif line[0] == "INS":
                obj.setdefault('insertions', [])
                obj['insertions'].append((line[1], line[2]))

    return obj


# def draw_center_of_mass(centers):
#     """
#     draws the center of mass of each unit

#     :param centers: list of coordinates of the centers of mass for each unit
#     :return:
#     """
#     for unit in centers:
#         pymol_center = [COLOR, 255, 0, 0, SPHERE] + list(centers[unit]) + [0.5]
#         cmd.load_cgo(pymol_center, 'center' + unit)
        
    
# def draw_distance_center_of_mass(centers):
#     """
#     draws the lines connecting the centers of mass 
#     (distances between two center of mass)
    
#     :param centers: list of coordinates of the centers of mass for each unit
#     :return:
#     """
#     for unit_1 in centers:
#         center_1 = centers[unit_1]
#         pymol.cmd.pseudoatom('center' + unit_1,
#                              pos=[center_1[0], center_1[1], center_1[2]],
#                              color="red", 
#                              name='center_1')
#         for unit_2 in centers:
#             if unit_1 != unit_2:
#                 center_2 = centers[unit_2]
#                 pymol.cmd.pseudoatom('center' + unit_2,
#                                      pos=[center_2[0], center_2[1], center_2[2]],
#                                      color="red", 
#                                      name='center_2')
#                 pymol.cmd.distance('center' + unit_1 + '////center_1',
#                                    'center' + unit_2 + '////center_2')

                            
# def draw_distance_center_mass_alpha(unit, center, alpha_c):
#     """
#     draws the lines between center of mass and alpha carbon in the unit 
#     (distance between center of mass and alpha carbons)

#     :param unit: string, name of the unit
#     :param center: coordinate representing the center of mass of the unit
#     :param alpha_c: list of coordinates representing the alpha carbon 
#                     of the unit
#     :return:
#     """
#     cmd.pseudoatom('center' + unit,
#                    pos=[center[0], center[1], center[2]],
#                    color="red", 
#                    name='center')
#     for i, c in enumerate(alpha_c):
#         coord = c.get_coord()
#         cmd.pseudoatom('c_alpha' + unit + str(i),
#                        pos=[coord[0], coord[1], coord[2]],
#                        color="red", 
#                        name='c_alpha')
#         cmd.distance('center' + unit + '////center',
#                      'c_alpha' + unit + str(i) + '////c_alpha') 


# def draw_vector(x, y):
#     """
#     draws a line (vector) between two point 

#     :param x: coordinates of the first point 
#     :param y: coordinates of the second point
#     :return:
#     """
#     id_ = str(uuid4())
#     cmd.pseudoatom('vx' + id_, 
#                    pos=x, 
#                    color="red", 
#                    name='x')
#     cmd.pseudoatom('vy' + id_, 
#                    pos=y, 
#                    color="red", 
#                    name='y')
#     cmd.distance('vx' + id_ + '////x',
#                  'vy' + id_ + '////y')