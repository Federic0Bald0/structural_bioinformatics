This project aims to describe a protein extracted from RepeatDB in its geometrical feature, providing also a way to visualize the results in:

    - JSON format
    - Visualization through PyMol 

The program is provided of a command line interface that allows the uset to interact with the tool using a series of yes or no questions. 
The functionality provided are:

    Center of mass, computes the centers of mass for each unit composing the protein. 

    Distance between center of mass, computes the distance between every center of mass of each unit in the protein. 

    Distance between center of mass and alpha carbon, computes the distance between the center of mass and each alpha carbon of the unit. 

    Handedness, for each unit is computed the handedness

    Twist, computes the rotation between consecutive units

    Curvature, it computes the horizontal angle between the vector connecting different units. 

    Pitch, it computes the vertical angle between the vector connecting different units. 

