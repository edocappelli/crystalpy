"""
Unittest for Photon Bunch class.
"""
import unittest


import numpy


from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon
from crystalpy.util.PhotonBunch import PhotonBunch

from numpy.testing import assert_array_almost_equal


class PhotonBunchTest(unittest.TestCase):
    def testFromList(self):


        npoint = 1000
        vx = numpy.zeros(npoint) + 0.0
        vy = numpy.zeros(npoint) + 1.0
        vz = numpy.zeros(npoint) + 0.0


        energy = numpy.zeros(npoint) + 3000.0

        photon_bunch1 = PhotonBunch()
        photon_bunch2 = PhotonBunch()

        photons_list = list()

        for i in range(npoint):

            photon = Photon(energy_in_ev=energy[i],
                                     direction_vector=Vector(vx[i],vy[i],vz[i]))

            photon_bunch1.addPhoton(photon)
            photons_list.append(photon)


        photon_bunch2.addPhotonsFromList(photons_list)

        energies = photon_bunch1.getArrayByKey("energies")
        for energy in energies:
            self.assertAlmostEqual(  energy, 3000.0)

        for i in range(len(photon_bunch1)):
            # print("checking photon %d "%i)
            self.assertTrue( photon_bunch1.getPhotonIndex(i) == photon_bunch2.getPhotonIndex(i) )



    def testChangePhotonValue(self):
        nphotons = 10


        bunch = PhotonBunch()
        list_of_photons = []
        for i in range(nphotons):
            photon = Photon(energy_in_ev=1000.0+i,
                                               direction_vector=Vector(1.0,0.0,0) )
            bunch.addPhoton(photon)
            list_of_photons.append(photon)


        for i in range(bunch.getNumberOfPhotons()):
            self.assertTrue( bunch.getPhotonIndex(i) == list_of_photons[i])

    def testToDictionary(self):
        self.assertTrue( isinstance( PhotonBunch().toDictionary(), dict))

    def testAddPhoton(self):

        photon1 = Photon(energy_in_ev=1000.0)
        photon2 = Photon(energy_in_ev=2000.0)

        bunch=PhotonBunch()

        self.assertTrue( bunch.getNumberOfPhotons() == 0)

        bunch.addPhoton(photon1)
        bunch.addPhoton(photon2)

        self.assertTrue( bunch.getNumberOfPhotons() == 2)

        bunch.addPhotonsFromList([photon1,photon2])

        self.assertTrue( bunch.getNumberOfPhotons() == 4)

        bunch.addBunch( bunch )

        self.assertTrue( bunch.getNumberOfPhotons() == 8)


    def testGetArrayByKey(self):

        photon = Photon(energy_in_ev=8000.0)
        bunch = PhotonBunch()

        for i in range(10):
            bunch.addPhoton(photon)

        assert_array_almost_equal( bunch.getArrayByKey("energies"), numpy.ones(10)*8000.0 )

        # for key in bunch.keys():
        #     print("array of ",key,bunch.getArrayByKey(key))

        # print( "<><>" , bunch.getArrayByKey("dd")  )



