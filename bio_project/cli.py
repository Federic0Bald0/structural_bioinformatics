# coding: utf-8

import compute_geometry

from PyInquirer import style_from_dict, Token, prompt
from PyInquirer import Validator, ValidationError

"""
TODO

beatify json print

review correctness questions

suppress bio python warning 

"""

style = style_from_dict({
    Token.QuestionMark: '#E91E63 bold',
    Token.Selected: '#673AB7 bold',
    Token.Instruction: '',  # default
    Token.Answer: '#2196f3 bold',
    Token.Question: '',
})

print ('Answer the question below to set the evaluation and visualization \
        paramenter for the protein')

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
        'message': 'Do you want to compute the distance between consecutive \
                    center of mass of the units?',
        'default': False
    },
    {
        'type': 'confirm',
        'name': 'c_alpha_distance',
        'message': 'Do you want to compute the distance betweem consecutive \
                    C-alpha in the units?',
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
        'message': 'Do you want to compute the twist (rotation) between \
                    consecutive units?',
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

]

answers = prompt(questions, style=style)

centers = compute_geometry.center_mass_unit()

if answers['center_of_mass']:
    if answers['json']:
        print centers 

if answers['center_of_mass_distance']:
    ris = compute_geometry.distance_center_of_mass(centers)
    if answers['json']: 
        print ris 

if answers['c_alpha_distance']:
    ris = compute_geometry.distance_alpha_c(centers)
    if answers['json']:
        print ris 

if answers['handedness']:
    ris = compute_geometry.handedness(centers)
    if answers['json']:
        print ris 

if answers['twist']:
    ris = compute_geometry.twist(centers)
    if answers['json']:
        print ris 

if answers['curvature'] or answers['pitch']:
    ris_c, ris_p = compute_geometry.curvature_pitch(centers)
    if answers['json'] and answers['curvature']:
        print ris_c 
    if answers['json'] and answers['pitch']:
        print ris_p 

