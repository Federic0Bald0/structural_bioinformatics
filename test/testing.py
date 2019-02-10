# coding: utf-8
import os

from unittest import TestCase

import bio_project.compute_geometry as compute_geometry
import bio_project.utils as utils

unit_test_dir = "test/test_unit/" 

v_1 = (5, 2, -3)
v_2 = (7, -8, -1)


class TestGeometry(TestCase):
    """
    class used to test the code
    """
    def test_distance(self):
        """
        test eclidean distance method
        :return:
        """ 
        ris = compute_geometry.distance(v_1, v_2)
        self.assertEquals(ris, 10.392304845413264)

    def test_angle(self):
        """
        test angle method
        :return:
        """ 
        ris = compute_geometry.angle(v_1, v_2)
        self.assertEquals(ris, 70.47273327273626)

    def test_atoi(self):
        """
        test atoi method
        :return:
        """ 
        ris = utils.atoi("42")
        # self.assertEquals(type(ris), int)
        self.assertEquals(ris, 42)

    def test_natural_key(self):
        """
        test natural_keys method
        :return:
        """ 
        ris = utils.natural_keys("int42")
        self.assertEquals(ris, ['int', 42])

    def test_center_of_mass(self):
        """
        test center_of_mass method
        :return:
        """ 
        ris = compute_geometry.center_mass_unit(False, unit_test_dir)
        self.assertEquals(ris['unit_27-48.pdb'], (-54.48145734501918,
                                                  30.31792712426065, 
                                                  15.1713042009396))
        self.assertEquals(ris['unit_49-72.pdb'], (-50.468922194033716,
                                                  34.10473405182134,
                                                  10.878357331957995))

    def test_distance_center_of_mass(self):
        """
        test distance_center_of_mass method
        :return:
        """ 
        centers = compute_geometry.center_mass_unit(False, unit_test_dir)
        ris = compute_geometry.distance_center_of_mass(centers,
                                                       False, 
                                                       unit_test_dir)
        self.assertEquals(ris['unit_27-48.pdb']['unit_49-72.pdb'], 
                          6.990689369755643)
        self.assertEquals(ris['unit_27-48.pdb']['unit_49-72.pdb'], 
                          ris['unit_49-72.pdb']['unit_27-48.pdb']) 

    def test_distance_alpha_c(self):
        """
        test distance_alpha_c method
        :return:
        """ 
        centers = compute_geometry.center_mass_unit(False, unit_test_dir)
        ris = compute_geometry.distance_alpha_c(centers,
                                                False,
                                                unit_test_dir)
        self.assertEqual(ris['unit_27-48.pdb']['CA1'], 8.609722419859281)

    def test_handedness(self):
        """
        test handedness method
        :return:
        """ 
        centers = compute_geometry.center_mass_unit(False, unit_test_dir)
        ris = compute_geometry.handedness(centers,
                                          False,
                                          unit_test_dir)
        self.assertEquals(ris['unit_27-48.pdb'], 'L')
        self.assertEquals(ris['unit_49-72.pdb'], 'L')

    def test_twist(self):
        """
        test twist method
        :return:
        """ 
        centers = compute_geometry.center_mass_unit(False, unit_test_dir)
        ris = compute_geometry.twist(centers,
                                     False,
                                     unit_test_dir)
        self.assertEquals(ris['unit_27-48.pdb_unit_49-72.pdb'],
                          0.9309720264364901)

    def test_pitch(self):
        """
        test pitch method
        :return:
        """ 
        centers = compute_geometry.center_mass_unit(False, unit_test_dir)
        ris = compute_geometry.pitch(centers,
                                     False,
                                     unit_test_dir)
        self.assertEquals(ris['unit_27-48.pdb_unit_49-72.pdb'],
                          11.498455699343008)

    def test_curvature(self):
        """
        test curvature method
        :return:
        """ 
        centers = compute_geometry.center_mass_unit(False, unit_test_dir)
        ris = compute_geometry.curvature(centers,
                                         False,
                                         unit_test_dir)
        self.assertEquals(ris['unit_27-48.pdb_unit_49-72.pdb'],
                          18.9001109414561)