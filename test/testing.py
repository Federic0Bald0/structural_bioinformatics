# coding: utf-8

from unittest import TestCase

import bio_project.compute_geometry as compute_geometry
import bio_project.utils as utils

unit_test_dir = "data/test_unit/"

v_1 = (5, 2, -3)
v_2 = (7, -8, -1)


class TestGeometry(TestCase):

    def test_distance(self): 
        ris = compute_geometry.distance(v_1, v_2)
        self.assertEquals(ris, 10.392304845413264)

    def test_angle(self):
        ris = compute_geometry.angle(v_1, v_2)
        self.assertEquals(ris, 70.47273327273626)

    def test_atoi(self):
        ris = utils.atoi("42")
        # self.assertEquals(type(ris), int)
        self.assertEquals(ris, 42)

    def test_natural_key(self):
        ris = utils.natural_keys("int42")
        self.assertEquals(ris, ['int', 42])

