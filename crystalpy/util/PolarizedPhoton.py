"""
This object contains a list of PolarizedPhoton objects, characterized by energy, direction vector and Stokes vector.
This object is used as input to and output from the passive crystal widget.
"""
from math import atan2

from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon
from crystalpy.util.StokesVector import StokesVector
from crystalpy.polarization.MuellerMatrix import MuellerMatrix

class PolarizedPhoton(Photon):
    """
    is a Photon object with a specified polarization state described by a Stokes vector.
    """
    def __init__(self, energy_in_ev, direction_vector, stokes_vector):
        """
        Constructor.
        :param energy_in_ev: photon energy in electron volts.
        :type energy_in_ev: float
        :param direction_vector: it doesn't have to be a unit vector, it gets normalized in Photon.__init__.
        :type direction_vector: crystalpy.util.Vector
        :param stokes_vector: Stokes vector describing the polarization state.
        :type stokes_vector: StokesVector
        """
        # self._energy_in_ev holds the energy in the base class.
        # self._unit_direction_vector holds the vector in the base class.
        self._stokes_vector = stokes_vector

        super(PolarizedPhoton, self).__init__(energy_in_ev, direction_vector)

    def set_unit_direction_vector(self, direction_vector):
        """
        :type direction_vector: Vector
        """
        self._unit_direction_vector = direction_vector.getNormalizedVector()

    def deviation(self):
        """
        the deviations are calculated supposing that the bunch moves along the y axis
        and considering a clockwise rotation as a positive deviation.
        """
        vector = self.unitDirectionVector().components()  # ndarray([x, y, z])
        deviation = atan2(vector[2], vector[1])

        return deviation

    def stokesVector(self):
        """
        :return: Stokes vector.
        """
        return self._stokes_vector

    def setStokesVector(self,stokes_vector):
        self._stokes_vector = stokes_vector

    def applyMuellerMatrix(self,mueller_matrix=MuellerMatrix()):
        s_in = self.stokesVector()
        s_out = mueller_matrix.calculate_stokes_vector( s_in )
        self.setStokesVector(s_out)


    def polarizationDegree(self):
        """
        :return: degree of circular polarization.
        """
        return self._stokes_vector.polarization_degree()

    def __eq__(self, candidate):
        """
        Determines if two polarized photons are identical (same energy, direction and polarization).
        :param candidate: Polarized photon to compare with.
        :return: True if equal otherwise False.
        """
        if ((self.energy() == candidate.energy() and
                self.unitDirectionVector() == candidate.unitDirectionVector()) and
                self.stokesVector() == candidate.stokesVector()):
            return True

        return False

    def __ne__(self, candidate):
        """
        Determines if two polarized photons are not identical (same energy, direction and polarization).
        :param candidate: Polarized photon to compare with.
        :return: True if not equal otherwise False.
        """
        return not (self == candidate)