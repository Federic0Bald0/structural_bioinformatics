# coding: utf-8
import os
import json
import pymol
import warnings

from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning  


import bio_project.compute_geometry as compute_geometry
from bio_project.utils import UnitSelection, parse_rdb

from PyInquirer import style_from_dict, Token, prompt
from PyInquirer import Validator, ValidationError


warnings.simplefilter('ignore', PDBConstructionWarning)

"""
TODO

beatify json print

review correctness questions

"""

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


style = style_from_dict({
    Token.QuestionMark: '#E91E63 bold',
    Token.Selected: '#673AB7 bold',
    Token.Instruction: '',  # default
    Token.Answer: '#2196f3 bold',
    Token.Question: '',
})

print('Answer the question below to set the evaluation and visualization'
      'paramenter for the protein:')

questions = [

    {
        'type': 'confirm',
        'name': 'center_of_mass',
        'message': 'Do you want to compute the centers of mass of the unit?',
        'default': False
    },
    {
        'type': 'confirm',
        'name': 'center_of_mass_distance',
        'message': ("Do you want to compute the distance between consecutive"
                    "center of mass of the units?"),
        'default': False
    },
    {
        'type': 'confirm',
        'name': 'c_alpha_distance',
        'message': ("Do you want to compute the distance betweem consecutive"
                    "C-alpha in the units?"),
        'default': False
    },
    {
        'type': 'confirm',
        'name': 'handedness',
        'message': 'Do you want to compute the handedness of unit?',
        'default': False
    },
    {
        'type': 'confirm',
        'name': 'twist',
        'message': ("Do you want to compute the twist (rotation) between"
                    "consecutive units?"),
        'default': False
    },
    {
        'type': 'confirm',
        'name': 'curvature',
        'message': 'Do you want to compute the cuvature?',
        'default': False
    },
    {
        'type': 'confirm',
        'name': 'pitch',
        'message': 'Do you want to compute the pitch?',
        'default': False
    },
    {   
        'type': 'confirm',
        'name': 'json',
        'message': 'Do you want to visualize the results in json format?',
        'default': False
    },
    {   
        'type': 'confirm',
        'name': 'draw',
        'message': 'Do you want to visualize the results with pymol?',
        'default': False
    },

]


answers = prompt(questions, style=style)

if answers['draw']:
    # pymol
    pymol.finish_launching()
    for unit in os.listdir(unit_dir):
        pymol.cmd.load(unit_dir + unit, unit)

centers = compute_geometry.center_mass_unit(answers['draw'], unit_dir)

if answers['center_of_mass']:
    if answers['json']:
        print(json.dumps(centers, indent=4, sort_keys=True))

if answers['center_of_mass_distance']:
    ris = compute_geometry.distance_center_of_mass(centers,
                                                   answers['draw'],
                                                   unit_dir)
    if answers['json']: 
        print(json.dumps(ris, indent=4, sort_keys=True)) 

if answers['c_alpha_distance']:
    ris = compute_geometry.distance_alpha_c(centers, answers['draw'], unit_dir)
    if answers['json']:
        print(json.dumps(ris, indent=4, sort_keys=True)) 

if answers['handedness']:
    ris = compute_geometry.handedness(centers, answers['draw'], unit_dir)
    if answers['json']:
        print(json.dumps(ris, indent=4, sort_keys=True)) 

if answers['twist']:
    ris = compute_geometry.twist(centers, answers['draw'], unit_dir)
    if answers['json']:
        print(json.dumps(ris, indent=4, sort_keys=True))

if answers['curvature']:
    ris = compute_geometry.curvature(centers, answers['draw'], unit_dir)
    if answers['json']:
        print(json.dumps(ris, indent=4, sort_keys=True))

if answers['pitch']: 
    ris = compute_geometry.pitch(centers, answers['draw'], unit_dir)
    if answers['json']:
        print(json.dumps(ris, indent=4, sort_keys=True)) 

