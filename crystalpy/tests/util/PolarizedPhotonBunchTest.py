"""
Unittest for Vector class.
"""
import unittest
import time

import numpy

from crystalpy.util.Vector import Vector
from crystalpy.util.StokesVector import StokesVector
from crystalpy.util.PolarizedPhotonBunch import PolarizedPhoton
from crystalpy.util.PolarizedPhotonBunch import PolarizedPhotonBunch


class PolarizedPhotonBunchTest(unittest.TestCase):
    def testFromArray(self):


        npoint = 1000
        vx = numpy.zeros(npoint) + 0.0
        vy = numpy.zeros(npoint) + 1.0
        vz = numpy.zeros(npoint) + 0.0

        s0 = numpy.zeros(npoint) + 1
        s1 = numpy.zeros(npoint) + 0
        s2 = numpy.zeros(npoint) + 1
        s3 = numpy.zeros(npoint) + 0

        energy = numpy.zeros(npoint) + 3000.0

        photon_bunch = PolarizedPhotonBunch([])

        t0 = time.time()
        photons_list = list()
        for i in range(npoint):

            photon = PolarizedPhoton(energy_in_ev=energy[i],
                                     direction_vector=Vector(vx[i],vy[i],vz[i]),
                                     stokes_vector=StokesVector([s0[i],s1[i],s2[i],s3[i]]))
            #photon_bunch.add(photon)
            photons_list.append(photon)


        photon_bunch.add(photons_list)
        print("Conversion in %f s"%(time.time() - t0))
        # print(photon_bunch.to_string())

        self.assertTrue(  1.0 == any(photon_bunch.get_array("s0"))  )
        self.assertTrue(  0.0 == any(photon_bunch.get_array("s1"))  )
        self.assertTrue(  1.0 == any(photon_bunch.get_array("s2"))  )
        self.assertTrue(  0.0 == any(photon_bunch.get_array("s3"))  )

        energies = photon_bunch.get_array("energies")
        for energy in energies:
            self.assertAlmostEqual(  energy, 3000.0)

        print(photon_bunch.photon_bunch)