"""
Unittest for GeometryType class.
"""
import unittest

from crystalpy.diffraction.GeometryType import BraggDiffraction, LaueDiffraction, \
                                      BraggTransmission, LaueTransmission
from crystalpy.diffraction.GeometryType import GeometryType


class GeometryTypeTest(unittest.TestCase):
    def testConstructor(self):
        geometry_type_description = "a geometry type"
        geometry_type = GeometryType(geometry_type_description)
        self.assertEqual(geometry_type.description(),
                         geometry_type_description)

    def testDescription(self):
        geometry_type_description = "a geometry type"
        geometry_type = GeometryType(geometry_type_description)
        self.assertEqual(geometry_type.description(),
                         geometry_type_description)

    def testEqualOperator(self):
        geometry_type_one = GeometryType("type one")
        geometry_type_two = GeometryType("type two")

        self.assertEqual(geometry_type_one,geometry_type_one)
        self.assertEqual(geometry_type_two,geometry_type_two)
        self.assertNotEqual(geometry_type_one,geometry_type_two)

    def testAllGeometryTypes(self):
        all_geometries = GeometryType.allGeometryTypes()

        self.assertIn(BraggDiffraction(), all_geometries)
        self.assertIn(LaueDiffraction(), all_geometries)
        self.assertIn(BraggTransmission(), all_geometries)
        self.assertIn(LaueTransmission(), all_geometries)


class BraggDiffractionTest(unittest.TestCase):
    def testConstructor(self):
        bragg_diffraction = BraggDiffraction()
        self.assertEqual(bragg_diffraction.description(),
                         "Bragg diffraction")


class LaueDiffractionTest(unittest.TestCase):
    def testConstructor(self):
        laue_diffraction = LaueDiffraction()
        self.assertEqual(laue_diffraction.description(),
                         "Laue diffraction")


class BraggTransmissionTest(unittest.TestCase):
    def testConstructor(self):
        bragg_transmission = BraggTransmission()
        self.assertEqual(bragg_transmission.description(),
                         "Bragg transmission")


class LaueTransmissionTest(unittest.TestCase):
    def testConstructor(self):
        laue_transmission = LaueTransmission()
        self.assertEqual(laue_transmission.description(),
                         "Laue transmission")