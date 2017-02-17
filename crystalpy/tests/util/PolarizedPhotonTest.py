"""
Unittest for PolarizedPhoton class.
"""

import unittest

from crystalpy.util.PolarizedPhoton import PolarizedPhoton
from crystalpy.util.Vector import Vector
from crystalpy.util.StokesVector import StokesVector

from numpy.testing import assert_array_almost_equal


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
        self.assertEqual(photon.energy(), 4000)

    def testDuplicate(self):
        ph1 = PolarizedPhoton(8000.0,Vector(2,4,5),StokesVector([1,2,3,4]))
        ph2 = ph1.duplicate()

        assert_array_almost_equal(ph1.stokesVector().components(),
                                  ph2.stokesVector().components() )

        assert_array_almost_equal(ph1.unitDirectionVector().components(),
                                  ph2.unitDirectionVector().components() )

    def testInheritatedMethods(self):
        ph1 = PolarizedPhoton(8000.0,Vector(2,4,5),StokesVector([1,2,3,4]))

        ph1.setUnitDirectionVector(Vector(1,0,0))
        ph1.setEnergy(9000)

        self.assertTrue( ph1.unitDirectionVector().components()[0] == 1)
        self.assertTrue( ph1.unitDirectionVector().components()[1] == 0)
        self.assertTrue( ph1.unitDirectionVector().components()[2] == 0)

        self.assertTrue( ph1.energy() == 9000.0)