This project aims to describe a protein extracted from RepeatDB in its geometrical feature, providing also a way to visualize the results in:

    - JSON format
    - Visualization through PyMol 

The program is provided of a command line interface that allows the uset to interact with the tool using a series of yes or no questions. 
The functionality provided are:

    Center of mass, computes the centers of mass for each unit composing the protein. The center of mass is computed as mean for each dimension, in respect 
    to the total mass of the unit, of the sum of the product of each atom mass for its component coordinate value.

    Distance between center of mass, computes the distance between every center of mass of each unit in the protein. This is achieved using eucledian distance between the set of coordinates representig the centers of mass.

    Distance between center of mass and alpha carbon, computes the distance between the center of mass and each alpha carbon of the unit. This is achieved using eucledian distance between the set of coordinates representig the center of mass and the alpha carbon.

    Handedness, for each unit is computed the handedness, that can be right-handed or left-handed. If the cross product of the vectores connecting the first and the second alpha carbon to the center of mass in the unit is positive the unit will be right-handed, otherwise it will be left-handed. 

    Twist, computes the rotation between consecutive units. This is achieved finding the two planes passing through the center of mass and the first alpha carbon and than computing the angle between this two.

    Curvature, it computes the horizontal angle between the vector connecting different units (through the center of mass). 

    Pitch

