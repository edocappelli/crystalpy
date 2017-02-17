"""
Unittest for DiffractionSetup class.
"""

import unittest

import numpy

from crystalpy.diffraction.DiffractionSetupSweeps import DiffractionSetupSweeps
from crystalpy.diffraction.GeometryType import BraggDiffraction
from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon


def diffractionSetup():


    diffraction_setup = DiffractionSetupSweeps(BraggDiffraction(),
                                         "Si",
                                         thickness=0.0001,
                                         miller_h=1,
                                         miller_k=1,
                                         miller_l=1,
                                         asymmetry_angle=10.0*numpy.pi/180,
                                         azimuthal_angle=0.0,
                                         energy_min=10000,
                                         energy_max=10000,
                                         energy_points=1,
                                         angle_deviation_min=-100.0e-6,
                                         angle_deviation_max=100e-6,
                                         angle_deviation_points=175)
    return diffraction_setup


class DiffractionSetupSweepsTest(unittest.TestCase):
    def testConstructor(self):
        diffraction_setup = diffractionSetup()
        self.assertIsInstance(diffraction_setup, DiffractionSetupSweeps)

        self.assertEqual(diffraction_setup._geometry_type,
                         BraggDiffraction())
        self.assertEqual(diffraction_setup._crystal_name,
                         "Si")
        self.assertEqual(diffraction_setup._thickness,
                         0.0001)
        self.assertEqual(diffraction_setup._miller_h,
                         1)
        self.assertEqual(diffraction_setup._miller_k,
                         1)
        self.assertEqual(diffraction_setup._miller_l,
                         1)
        self.assertEqual(diffraction_setup._asymmetry_angle,
                         10*numpy.pi/180)
        self.assertEqual(diffraction_setup._azimuthal_angle,
                         0)
        self.assertEqual(diffraction_setup.energyMin(),
                         10000)
        self.assertEqual(diffraction_setup.energyMax(),
                         10000)
        self.assertEqual(diffraction_setup.energyPoints(),
                         1)

        self.assertAlmostEqual(diffraction_setup.angleDeviationMin(),
                               -100.0e-6)
        self.assertAlmostEqual(diffraction_setup.angleDeviationMax(),
                               100.0e-6)
        self.assertEqual(diffraction_setup.angleDeviationPoints(),
                         175)

    def testGeometryType(self):
        diffraction_setup = diffractionSetup()
        self.assertEqual(diffraction_setup.geometryType(),
                         BraggDiffraction())

    def testCrystalName(self):
        diffraction_setup = diffractionSetup()
        self.assertEqual(diffraction_setup.crystalName(),
                         "Si")

    def testThickness(self):
        diffraction_setup = diffractionSetup()
        self.assertEqual(diffraction_setup.thickness(),
                         0.0001)

    def testMillerH(self):
        diffraction_setup = diffractionSetup()
        self.assertEqual(diffraction_setup.millerH(),
                         1)

    def testMillerK(self):
        diffraction_setup = diffractionSetup()
        self.assertEqual(diffraction_setup.millerK(),
                         1)

    def testMillerL(self):
        diffraction_setup = diffractionSetup()
        self.assertEqual(diffraction_setup.millerL(),
                         1)

    def testAsymmetryAngle(self):
        diffraction_setup = diffractionSetup()
        self.assertEqual(diffraction_setup.asymmetryAngle(),
                         10*numpy.pi/180)

    def testAzimuthalAngle(self):
        diffraction_setup = diffractionSetup()
        self.assertEqual(diffraction_setup.azimuthalAngle(),
                         0)

    def testEnergyMin(self):
        diffraction_setup = diffractionSetup()
        self.assertEqual(diffraction_setup.energyMin(),
                         10000)

    def testEnergyMax(self):
        diffraction_setup = diffractionSetup()
        self.assertEqual(diffraction_setup.energyMax(),
                         10000)

    def testEnergyMin(self):
        diffraction_setup = diffractionSetup()
        self.assertEqual(diffraction_setup.energyMin(),
                         10000)

    def testAngleDeviationMin(self):
        diffraction_setup = diffractionSetup()
        self.assertAlmostEqual(diffraction_setup.angleDeviationMin(),
                               -100.0e-6)

    def testAngleDeviationMax(self):
        diffraction_setup = diffractionSetup()
        self.assertAlmostEqual(diffraction_setup.angleDeviationMax(),
                               100.0e-6)

    def testAngleDeviationPoints(self):
        diffraction_setup = diffractionSetup()
        self.assertEqual(diffraction_setup.angleDeviationPoints(),
                         175)

    def testUnitcellVolume(self):
        diffraction = diffractionSetup()

        unitcell_volume = diffraction.unitcellVolume()
        self.assertAlmostEqual(unitcell_volume, 160.1649322509)

    def testAsInfoDictionary(self):
        diffraction_setup = diffractionSetup()

        info_dict = diffraction_setup.toDictionary()

        self.assertEqual(info_dict["Geometry Type"],
                         "Bragg diffraction")
        self.assertEqual(info_dict["Crystal Name"],
                         "Si")
        self.assertEqual(info_dict["Thickness"],
                         "0.0001")
        self.assertEqual(info_dict["Miller indices (h,k,l)"],
                         "(1,1,1)")
        self.assertEqual(info_dict["Asymmetry Angle"],
                         str(10*numpy.pi/180))
        self.assertEqual(info_dict["Azimuthal Angle"],
                         str(0.0))
        self.assertEqual(info_dict["Minimum energy"],
                         "10000.0")
        self.assertEqual(info_dict["Maximum energy"],
                         "10000.0")
        self.assertEqual(info_dict["Number of energy points"],
                         "1")
        self.assertEqual(info_dict["Angle deviation minimum"],
                         "-1.00e-04")
        self.assertEqual(info_dict["Angle deviation maximum"],
                         "1.00e-04")
        self.assertEqual(info_dict["Angle deviation points"],
                         "175")

    def testOperatorEqual(self):
        diffraction_setup_one = diffractionSetup()
        diffraction_setup_two = DiffractionSetupSweeps(BraggDiffraction(),
                                                 "Diamond",
                                                 thickness=0.001,
                                                 miller_h=1,
                                                 miller_k=1,
                                                 miller_l=1,
                                                 asymmetry_angle=10*numpy.pi/180,
                                                 azimuthal_angle=0,
                                                 energy_min=8000,
                                                 energy_max=8000,
                                                 energy_points=1,
                                                 angle_deviation_min= -100.0e-6,
                                                 angle_deviation_max=100e-6,
                                                 angle_deviation_points=175)


        self.assertTrue(diffraction_setup_one == diffractionSetup())
        self.assertFalse(diffraction_setup_one == diffraction_setup_two)


    def testOperatorNotEqual(self):
        diffraction_setup_one = diffractionSetup()
        diffraction_setup_two = DiffractionSetupSweeps(BraggDiffraction(),
                                                 "Diamond",
                                                 thickness=0.001,
                                                 miller_h=1,
                                                 miller_k=1,
                                                 miller_l=1,
                                                 asymmetry_angle=11,
                                                 azimuthal_angle=0.0,
                                                 energy_min=8000,
                                                 energy_max=8000,
                                                 energy_points=1,
                                                 angle_deviation_min= -100.0e-6,
                                                 angle_deviation_max=100e-6,
                                                 angle_deviation_points=175)


        self.assertTrue(diffraction_setup_one != diffraction_setup_two)
        self.assertFalse(diffraction_setup_one != diffractionSetup())

    def testClone(self):
        diffraction_setup = diffractionSetup()
        clone = diffraction_setup.clone()

        self.assertEqual(diffraction_setup, clone)
        self.assertIsNot(diffraction_setup, clone)

    #
    # these are duplicated in DiffractionSetup, but do not harm...
    #
    def testAngleDeviationGrid(self):
        diffraction_setup = diffractionSetup()
        self.assertAlmostEqual(numpy.linalg.norm(diffraction_setup.angleDeviationGrid()-numpy.linspace(-100.0e-6, 100.0e-6, 175)),
                               0.0)

    def testAngleBragg(self):
        diffraction = diffractionSetup()
        angle_bragg = diffraction.angleBragg(energy=8000)
        self.assertAlmostEqual(angle_bragg, 0.249732328921)

    def testF0(self):
        diffraction = diffractionSetup()

        f_0 = diffraction.F0(energy=8000)
        self.assertAlmostEqual(f_0, 114.08416+2.7188j)

    def testFH(self):
        diffraction = diffractionSetup()

        f_h = diffraction.FH(energy=8000)
        self.assertAlmostEqual(f_h, 44.54356349760925-41.82476349760927j)

    def testFH_bar(self):
        diffraction = diffractionSetup()

        f_h_bar = diffraction.FH_bar(energy=8000)
        self.assertAlmostEqual(f_h_bar, 41.82476349760923+44.54356349760926j)

    def testDSpacing(self):
        diffraction = diffractionSetup()

        d_spacing = diffraction.dSpacing()
        self.assertAlmostEqual(d_spacing, 3.135416288633)

    def testNormalBragg(self):
        diffraction = diffractionSetup()
        bragg_normal = diffraction.normalBragg()
        self.assertAlmostEqual(bragg_normal.getNormalizedVector().components()[0] , 0.0)
        self.assertAlmostEqual(bragg_normal.getNormalizedVector().components()[1] , numpy.sin(10*numpy.pi/180))
        self.assertAlmostEqual(bragg_normal.getNormalizedVector().components()[2] , numpy.cos(10*numpy.pi/180))

    def testNormalSurface(self):
        diffraction = diffractionSetup()

        surface_normal = diffraction.normalSurface()
        self.assertAlmostEqual(surface_normal.getNormalizedVector().components()[0] , 0.0)
        self.assertAlmostEqual(surface_normal.getNormalizedVector().components()[1] , 0.0)
        self.assertAlmostEqual(surface_normal.getNormalizedVector().components()[2] , 1.0)

    def testIncomingPhotonDirection(self):
        diffraction = diffractionSetup()

        photon_direction = diffraction.incomingPhotonDirection(8000, 0.01)
        self.assertEqual(photon_direction,
                         Vector(0.,0.90717943,-0.42074397))

    def testDeviationOfIncomingPhoton(self):
        diffraction = diffractionSetup()

        for energy in [2500, 6000, 8000, 15000, 22000, 30000]:
            for test_deviation in [0.01, 0.03, 0.5, -0.1, -0.9, 0.00001, -0.0007]:
                photon_direction = diffraction.incomingPhotonDirection(energy, test_deviation)
                photon = Photon(energy, photon_direction)
                deviation = diffraction.deviationOfIncomingPhoton(photon)
                self.assertAlmostEqual(test_deviation, deviation)
