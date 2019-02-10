# coding: utf-8
import os
import json
import pymol
import warnings

from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning  


import bio_project.compute_geometry as compute_geometry
from bio_project.utils import UnitSelection, parse_rdb, natural_keys

from PyInquirer import style_from_dict, Token, prompt
from PyInquirer import Validator, ValidationError


warnings.simplefilter('ignore', PDBConstructionWarning)

pdb_f = "data/2z7x.pdb"
rdb_f = "data/rdb/2z7xB.db"
unit_dir = "data/unit/"
unit_dir_list = os.listdir(unit_dir)
unit_dir_list.sort(key=natural_keys) 
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
        'message': 'Compute the centers of mass of the unit?',
        'default': False
    },
    {
        'type': 'confirm',
        'name': 'center_of_mass_distance',
        'message': ("Compute the distance between centers of mass?"),
        'default': False
    },
    {
        'type': 'confirm',
        'name': 'c_alpha_distance',
        'message': ("Compute the distance betweem center of mass"
                    " and alpha carbon in the units?"),
        'default': False
    },
    {
        'type': 'confirm',
        'name': 'handedness',
        'message': 'Compute the handedness?',
        'default': False
    },
    {
        'type': 'confirm',
        'name': 'twist',
        'message': ("Compute the twist (rotation) between consecutive units?"),
        'default': False
    },
    {
        'type': 'confirm',
        'name': 'curvature',
        'message': 'Compute the cuvature?',
        'default': False
    },
    {
        'type': 'confirm',
        'name': 'pitch',
        'message': 'Compute the pitch?',
        'default': False
    },
    {   
        'type': 'confirm',
        'name': 'json',
        'message': 'Display the results in json format?',
        'default': False
    },
    {   
        'type': 'confirm',
        'name': 'draw',
        'message': 'Display the results with pymol?',
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
    else:
        for unit in unit_dir_list:
            keys = natural_keys(unit)
            print("Center of mass of " + unit.split('.')[0] + " is: " +
                  str(centers[unit]))

if answers['center_of_mass_distance']:
    ris = compute_geometry.distance_center_of_mass(centers,
                                                   answers['draw'],
                                                   unit_dir)
    if answers['json']: 
        print(json.dumps(ris, indent=4, sort_keys=True)) 
    else:
        for unit_1 in unit_dir_list:
            for unit_2 in unit_dir_list:
                if unit_1 != unit_2:
                    print("The distance between the center of mass of" + 
                          unit_1.split('.')[0] + " and " +
                          unit_2.split('.')[0] + " is: " + 
                          str(ris[unit_1][unit_2])) 

if answers['c_alpha_distance']:
    ris = compute_geometry.distance_alpha_c(centers, answers['draw'], unit_dir)
    if answers['json']:
        print(json.dumps(ris, indent=4, sort_keys=True)) 
    else:
        for unit in unit_dir_list:
            for ca in ris[unit]:
                dist = str(ris[unit][ca])
                print("The distance between " + unit.split('.')[0] + 
                      " and alpha carbon " + str(natural_keys(ca)[1]) + 
                      " is: " + dist)

if answers['handedness']:
    ris = compute_geometry.handedness(centers, answers['draw'], unit_dir)
    if answers['json']:
        print(json.dumps(ris, indent=4, sort_keys=True), unit_dir) 
    else:
        for unit in unit_dir_list:
            keys = natural_keys(unit)
            hand = 'left-handed' if ris[unit] == 'L' else 'right-handed'
            print("The " + unit.split('.')[0] + " is " + hand)

if answers['twist']:
    ris = compute_geometry.twist(centers, answers['draw'], unit_dir)
    if answers['json']:
        print(json.dumps(ris, indent=4, sort_keys=True))
    else:
        for i in range(len(unit_dir_list)-1):
            unit_1 = unit_dir_list[i]
            unit_2 = unit_dir_list[i+1]
            twist = str(ris[unit_1 + '_' + unit_2])
            print("The twist between " + 
                  unit_1.split('.')[0] + " and " +
                  unit_2.split('.')[0] + " is: " + 
                  twist)

if answers['curvature']:
    ris = compute_geometry.curvature(centers, answers['draw'], unit_dir)
    if answers['json']:
        print(json.dumps(ris, indent=4, sort_keys=True))
    else:
        for i in range(len(unit_dir_list)-2):
            unit_1 = unit_dir_list[i]
            unit_2 = unit_dir_list[i+1]
            curv = str(ris[unit_1 + '_' + unit_2])
            print("The curvature between " + 
                  unit_1.split('.')[0] + " and " +
                  unit_2.split('.')[0] + " is: " + 
                  curv)

if answers['pitch']: 
    ris = compute_geometry.pitch(centers, answers['draw'], unit_dir)
    if answers['json']:
        print(json.dumps(ris, indent=4, sort_keys=True)) 
    else:
        for i in range(len(unit_dir_list)-2):
            unit_1 = unit_dir_list[i]
            unit_2 = unit_dir_list[i+1]
            pitch = str(ris[unit_1 + '_' + unit_2])
            print("The pitch between " + 
                  unit_1.split('.')[0] + " and " +
                  unit_2.split('.')[0] + " is: " + 
                  pitch)