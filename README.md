This project aims to describe a protein extracted from RepeatDB in its geometrical feature, providing also a way to visualize the results in:

    - JSON format
    - PyMol 

The program is provided of a command line interface that allows the user to interact with through a series of yes or no questions. 
The functionality provided are:

    Center of mass, computes the centers of mass for each unit composing the protein. 

    Distance between centers of mass, computes the distance between every center of mass of each unit in the protein. 

    Distance between center of mass and alpha carbons, computes the distance between the center of mass and each alpha carbon of the unit. 

    Handedness, for each unit is computed the handedness, every unit can be left-handed or right-handed.

    Twist, computes the rotation between consecutive units.

    Curvature, computes the horizontal angle between the vectors connecting consecutive units. 

    Pitch, computes the vertical angle between the vectors connecting consecutive units. 

