"""
Unittest for Vector class.
"""
import unittest
import time

import numpy

from numpy.testing import assert_almost_equal

from crystalpy.util.Vector import Vector
from crystalpy.util.StokesVector import StokesVector
from crystalpy.util.PolarizedPhoton import PolarizedPhoton
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

        photons_list = list()
        for i in range(npoint):

            photon = PolarizedPhoton(energy_in_ev=energy[i],
                                     direction_vector=Vector(vx[i],vy[i],vz[i]),
                                     stokes_vector=StokesVector([s0[i],s1[i],s2[i],s3[i]]))
            photons_list.append(photon)


        photon_bunch.addPhotonsFromList(photons_list)

        # print("<><><>",photon_bunch.toString())
        # print("<><><>",photon_bunch.toDictionary())

        self.assertTrue(  1.0 == any(photon_bunch.getArrayByKey("s0"))  )
        self.assertTrue(  0.0 == any(photon_bunch.getArrayByKey("s1"))  )
        self.assertTrue(  1.0 == any(photon_bunch.getArrayByKey("s2"))  )
        self.assertTrue(  0.0 == any(photon_bunch.getArrayByKey("s3"))  )

        energies = photon_bunch.getArrayByKey("energies")
        for energy in energies:
            self.assertAlmostEqual(  energy, 3000.0)



    def testChangePhotonValue(self):
        nphotons = 10

        from crystalpy.util.Vector import Vector
        from crystalpy.util.StokesVector import StokesVector

        bunch = PolarizedPhotonBunch([])
        for i in range(nphotons):
            polarized_photon = PolarizedPhoton(energy_in_ev=1000.0+i,
                                               direction_vector=Vector(0,1.0,0),
                                               stokes_vector=StokesVector([1.0,0,1.0,0]))
            bunch.addPhoton(polarized_photon)

        # photon5_stokes = bunch.get_photon_index(5).stokesVector().get_array(numpy=True)
        # print("photon 5 stokes ",photon5_stokes)

        photon5 = bunch.getPhotonIndex(5)

        photon5.setStokesVector(StokesVector([1,0,0,0]))
        photon5_stokes_new = bunch.getPhotonIndex(5).stokesVector().components()


        # print("photon 5 stokes new ",photon5_stokes_new)

        assert_almost_equal(photon5_stokes_new,numpy.array([1.0,0,0,0]))




