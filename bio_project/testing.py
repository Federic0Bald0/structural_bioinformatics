# coding: utf-8

from unittest import TestCase

import bio_project.compute_geometry as compute_geometry

unit_test_dir = "data/test_unit/"


class TestGeometry(TestCase):

    def test_center_of_mass(self): 
        ris = compute_geometry.center_mass_unit(False, unit_test_dir)
        print(ris)
        self.assertEquals(ris, 3)

    def test_distance_center_of_mass(self):
        pass

    def test_distance_alpha_c(self):
        pass

    def test_handedness(self):
        pass

    def test_twist(self):
        pass

    def test_curvature(self):
        pass

    def test_pitch(self):
        pass

ts = TestGeometry()
ts.test_center_of_mass()