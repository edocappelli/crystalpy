"""
Unittest for DiffractionSetup class.
"""

import unittest

import numpy

from crystalpy.diffraction.DiffractionSetup import DiffractionSetup
from crystalpy.diffraction.GeometryType import BraggDiffraction
from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon


def diffractionSetup():
    directions = [Vector(0,0,1.0/float(n)) for n in range(1,176)]
    photons = [Photon(10000, direction) for direction in directions]
    diffraction_setup = DiffractionSetup(BraggDiffraction(),
                                         "Si",
                                         thickness=0.0001,
                                         miller_h=1,
                                         miller_k=1,
                                         miller_l=1,
                                         asymmetry_angle=10*numpy.pi/180,
                                         azimuthal_angle=0.0,)
                                         # incoming_photons=photons)
    return diffraction_setup


class DiffractionSetupTest(unittest.TestCase):
    def testConstructor(self):
        diffraction_setup = diffractionSetup()
        self.assertIsInstance(diffraction_setup, DiffractionSetup)

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
    def testAzymuthalAngle(self):
        diffraction_setup = diffractionSetup()
        self.assertEqual(diffraction_setup.azimuthalAngle(),
                         0)

    # def testEnergyMin(self):
    #     diffraction_setup = diffractionSetup()
    #     self.assertEqual(diffraction_setup.energyMin(),
    #                      10000)
    #
    # def testEnergyMax(self):
    #     diffraction_setup = diffractionSetup()
    #     self.assertEqual(diffraction_setup.energyMax(),
    #                      10000)
    #
    # def testEnergyMin(self):
    #     diffraction_setup = diffractionSetup()
    #     self.assertEqual(diffraction_setup.energyMin(),
    #                      10000)

    # @unittest.expectedFailure
    # def testAngleDeviationMin(self):
    #     diffraction_setup = diffractionSetup()
    #     self.assertAlmostEqual(diffraction_setup.angleDeviationMin(),
    #                            -100.0e-6)
    # @unittest.expectedFailure
    # def testAngleDeviationMax(self):
    #     diffraction_setup = diffractionSetup()
    #     self.assertAlmostEqual(diffraction_setup.angleDeviationMax(),
    #                            100.0e-6)

    # def testAngleDeviationPoints(self):
    #     diffraction_setup = diffractionSetup()
    #     self.assertEqual(diffraction_setup.angleDeviationPoints(),
    #                      175)

    # @unittest.expectedFailure
    # def testAngleDeviationGrid(self):
    #     diffraction_setup = diffractionSetup()
    #     self.assertAlmostEqual(numpy.linalg.norm(diffraction_setup.angleDeviationGrid()-numpy.linspace(-100.0e-6, 100.0e-6, 175)),
    #                            0.0)

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

        bragg_normal = diffraction.normalBragg(return_normalized=True)
        # print("<><>",bragg_normal.components())
        # [ 0.          0.17364818  0.98480775]

        self.assertAlmostEqual (bragg_normal.components()[0],0.)
        self.assertAlmostEqual (bragg_normal.components()[1],0.17364818)
        self.assertAlmostEqual (bragg_normal.components()[2],0.98480775)

    def testNormalSurface(self):
        diffraction = diffractionSetup()

        surface_normal = diffraction.normalSurface()
        self.assertEqual(surface_normal,
                         Vector(0,0,1))

    def testIncomingPhotonDirection(self):
        diffraction = diffractionSetup()

        photon_direction = diffraction.incomingPhotonDirection(8000, 0)
        self.assertAlmostEqual (photon_direction.components()[0],0.)
        self.assertAlmostEqual (photon_direction.components()[1],0.91134144)
        self.assertAlmostEqual (photon_direction.components()[2],-0.41165129)



        photon_direction = diffraction.incomingPhotonDirection(8000, 0.01)
        self.assertAlmostEqual (photon_direction.components()[0],0.)
        self.assertAlmostEqual (photon_direction.components()[1],0.90717943)
        self.assertAlmostEqual (photon_direction.components()[2],-0.42074397)


    # def testDeviationOfIncomingPhoton(self):
    #     diffraction = diffractionSetup()
    #
    #     for energy in [2500, 6000, 8000, 15000, 22000, 30000]:
    #         for test_deviation in [0.01, 0.03, 0.5, -0.1, -0.9, 0.00001, -0.0007]:
    #             photon_direction = diffraction.incomingPhotonDirection(energy, test_deviation)
    #             photon = Photon(energy, photon_direction)
    #             deviation = diffraction.deviationOfIncomingPhoton(photon)
    #             self.assertAlmostEqual(test_deviation, deviation)

    def testUnitcellVolume(self):
        diffraction = diffractionSetup()

        unitcell_volume = diffraction.unitcellVolume()
        self.assertAlmostEqual(unitcell_volume, 160.1649322509)

    # @unittest.expectedFailure
    def testAsInfoDictionary(self):
        diffraction_setup = diffractionSetup()

        info_dict = diffraction_setup.asInfoDictionary()

        self.assertEqual(info_dict["Geometry Type"],
                         "Bragg diffraction")
        self.assertEqual(info_dict["Crystal Name"],
                         "Si")
        self.assertEqual(info_dict["Thickness"],
                         "0.0001")
        self.assertEqual(info_dict["Miller indices (h,k,l)"],
                         "(1,1,1)")
        self.assertEqual(info_dict["Asymmetry Angle"],
                         str(10.0*numpy.pi/180))
        self.assertEqual(info_dict["Azimuthal Angle"],
                         str(0.0))
        # self.assertEqual(info_dict["Minimum energy"],
        #                  "10000.0")
        # self.assertEqual(info_dict["Maximum energy"],
        #                  "10000.0")
        # self.assertEqual(info_dict["Number of energy points"],
        #                  "1")
        # self.assertEqual(info_dict["Angle deviation minimum"],
        #                  "-1.00e-04")
        # self.assertEqual(info_dict["Angle deviation maximum"],
        #                  "1.00e-04")
        # self.assertEqual(info_dict["Angle deviation points"],
        #                  "175")

    # @unittest.expectedFailure
    def testOperatorEqual(self):
        diffraction_setup_one = diffractionSetup()
        diffraction_setup_two = DiffractionSetup(BraggDiffraction(),
                                                 "Diamond",
                                                 thickness=0.001,
                                                 miller_h=1,
                                                 miller_k=1,
                                                 miller_l=1,
                                                 asymmetry_angle=11,
                                                 azimuthal_angle=0.0,)
                                                 # incoming_photons=None)


        self.assertTrue(diffraction_setup_one == diffractionSetup())
        self.assertFalse(diffraction_setup_one == diffraction_setup_two)

    # @unittest.expectedFailure
    def testOperatorNotEqual(self):
        diffraction_setup_one = diffractionSetup()
        diffraction_setup_two = DiffractionSetup(BraggDiffraction(),
                                                 "Diamond",
                                                 thickness=0.001,
                                                 miller_h=1,
                                                 miller_k=1,
                                                 miller_l=1,
                                                 asymmetry_angle=10*numpy.pi/180,
                                                 azimuthal_angle=0.0,)
                                                 # incoming_photons=None)


        self.assertTrue(diffraction_setup_one != diffraction_setup_two)
        self.assertFalse(diffraction_setup_one != diffractionSetup())

    def testClone(self):
        diffraction_setup = diffractionSetup()
        clone = diffraction_setup.clone()

        self.assertEqual(diffraction_setup, clone)
        self.assertIsNot(diffraction_setup, clone)


    def testGeometry(self):
        diffraction_setup = diffractionSetup()
        self.assertTrue( diffraction_setup.asInfoDictionary()["Geometry Type"] == "Bragg diffraction")

