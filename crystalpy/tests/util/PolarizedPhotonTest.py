"""
Unittest for PolarizedPhoton class.
"""

import unittest

from crystalpy.util.PolarizedPhoton import PolarizedPhoton
from crystalpy.util.Vector import Vector
from crystalpy.util.StokesVector import StokesVector


class PolarizedPhotonTest(unittest.TestCase):
    def testConstructor(self):
        photon = PolarizedPhoton(energy_in_ev=8000,
                                 direction_vector=Vector(0.0,1.0,0.0),
                                 stokes_vector=StokesVector( [1.0,0.0,1.0,0.0] ))

        self.assertIsInstance(photon, PolarizedPhoton)
        self.assertTrue(photon.unitDirectionVector() == Vector(0.0,1.0,0.0))
        self.assertTrue(photon.stokesVector() == StokesVector( [1.0,0.0,1.0,0.0] ))

    def testEnergy(self):
        photon = PolarizedPhoton(4000, Vector(0, 0, 1), StokesVector([1.0,0.0,1.0,0.0] ))
        self.assertEqual(photon.energy(), 8000)
