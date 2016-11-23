"""
Unittest for Photon class.
"""

import unittest

from orangecontrib.crystal.util.Photon import Photon
from orangecontrib.crystal.util.Vector import Vector

class PhotonTest(unittest.TestCase):
    def testConstructor(self):
        photon = Photon(4000, Vector(0, 0, 1))

        self.assertIsInstance(photon, Photon)
        self.assertEqual(photon.energy(), 4000)
        self.assertTrue(photon.unitDirectionVector() == Vector(0, 0, 1))

    def testEnergy(self):
        photon = Photon(4000, Vector(0, 0, 1))
        self.assertEqual(photon.energy(), 4000)

    def testWavelength(self):
        # Test data in eV : m.
        test_data = {   3 : 413.28 * 10e-9,
                        4 : 309.96 * 10e-9,
                        8 : 154.98 * 10e-9,
                     5000 : 2.4797 * 10e-10,
                    10000 : 1.2398 * 10e-10}

        for energy, wavelength in test_data.items():
            photon = Photon(energy, Vector(0, 0, 1))
            self.assertAlmostEqual(photon.wavelength(),
                                   wavelength, 2)

    def testWavenumber(self):
        # Test data in eV : m^-1.
        test_data = {   3 : 15203192.8,
                        4 : 20270923.76,
                        8 : 40541847.5,
                     5000 : 25338654707.5,
                    10000 : 50677309415}


        for energy, wavenumber in test_data.items():
            photon = Photon(energy, Vector(0, 0, 1))
            self.assertAlmostEqual(photon.wavenumber(),
                                   wavenumber, 1)

    def testWavevector(self):
        direction = Vector(0, 0, 1)
        photon = Photon(5000.0, direction)

        wavevector = photon.wavevector()

        self.assertAlmostEqual(wavevector.norm(),
                               25338654707.5, 1)

        self.assertEqual(wavevector.getNormalizedVector(),
                         direction)

    def testUnitDirectionVector(self):
        photon = Photon(4000, Vector(0, 0, 5))

        self.assertTrue(photon.unitDirectionVector() == Vector(0, 0, 1))

    def testOperatorEqual(self):
        photon_one = Photon(4000, Vector(0, 0, 5))
        photon_two = Photon(4000, Vector(0, 1, 1))
        photon_three = Photon(2000, Vector(0, 0, 5))

        self.assertTrue(photon_one == photon_one)
        self.assertFalse(photon_one == photon_two)
        self.assertFalse(photon_one == photon_three)
        self.assertFalse(photon_two == photon_three)

    def testOperatorNotEqual(self):
        photon_one = Photon(4000, Vector(0, 0, 5))
        photon_two = Photon(4000, Vector(0, 1, 1))
        photon_three = Photon(2000, Vector(0, 0, 5))

        self.assertFalse(photon_one != photon_one)
        self.assertTrue(photon_one != photon_two)
        self.assertTrue(photon_one != photon_three)
        self.assertTrue(photon_two != photon_three)