"""
Represents a photon.
Except for energy all units are in SI. Energy is in eV.
"""

from crystalpy.util.Vector import Vector
import scipy.constants as codata
import numpy


class Photon(object):

    def __init__(self, energy_in_ev=1000.0, direction_vector=Vector(0.0,1.0,0.0)):
        """
        Constructor.
        :param energy_in_ev: Photon energy in eV.
        :param direction_vector: The direction of the photon as Vector.
        """
        self._energy_in_ev = float(energy_in_ev)
        self._unit_direction_vector = direction_vector.getNormalizedVector()



    def duplicate(self):
        return Photon( self.energy(), self.unitDirectionVector() )

    def energy(self):
        """
        :return: Energy in eV.
        """
        return self._energy_in_ev

    def setEnergy(self,value):
        self._energy_in_ev = value

    def wavelength(self):
        """
        :return: The photon wavelength in meter.
        """
        E_in_Joule = self.energy() * codata.e # elementary_charge
        # Wavelength in meter
        wavelength = (codata.c * codata.h / E_in_Joule)
        return wavelength

    def wavenumber(self):
        """
        :return: Wavenumber in m^-1.
        """
        return (2.0 * numpy.pi) / self.wavelength()

    def wavevector(self):
        """
        :return: Photon wavevector in m^-1.
        """
        return self.unitDirectionVector().scalarMultiplication(self.wavenumber())

    def unitDirectionVector(self):
        """
        :return: Photon direction.
        """
        return self._unit_direction_vector

    def setUnitDirectionVector(self,vector=Vector(0,1,0)):
        """
        :return: Photon direction.
        """
        self._unit_direction_vector = vector.getNormalizedVector()


    def deviation(self):
        """
        the deviations are calculated supposing that the bunch moves along the y axis
        """
        vector = self.unitDirectionVector().components()  # ndarray([x, y, z])
        deviation = numpy.arctan2(vector[2], vector[1])

        return deviation

    def __eq__(self, candidate):
        """
        Determines if two photons are identical (same energy and direction).
        :param candidate: Photon to compare with.
        :return: True if equal otherwise False.
        """
        if (self.energy() == candidate.energy() and
                self.unitDirectionVector() == candidate.unitDirectionVector()):
            return True

        return False

    def __ne__(self, candidate):
        """
        Determines if two photons are not identical (same energy and direction).
        :param candidate: Photon to compare with.
        :return: True if not equal otherwise False.
        """
        return not (self == candidate)
