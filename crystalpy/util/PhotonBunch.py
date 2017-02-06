"""
This object contains a list of PolarizedPhoton objects, characterized by energy, direction vector and Stokes vector.
This object is used as input to and output from the passive crystal widget.
"""
from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon
from crystalpy.polarization.StokesVector import StokesVector
import numpy as np
from math import atan2


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
        deviation = atan2(-vector[2], vector[1])

        return deviation

    def stokesVector(self):
        """
        :return: Stokes vector.
        """
        return self._stokes_vector

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


class PhotonBunch(object):
    """
    is a collection of PolarizedPhoton objects, making up the photon beam.
    """
    def __init__(self, polarized_photons):
        """
        :param polarized_photons: bunch of PolarizedPhoton objects.
        :type polarized_photons: list(PolarizedPhoton, PolarizedPhoton, ...)
        """
        self.photon_bunch = polarized_photons
        self._set_dict()

    def add(self, to_be_added):
        """
        :param to_be_added: PolarizedPhoton object(s) to be added to the bunch.
        :type to_be_added: PolarizedPhoton or list(PolarizedPhoton) of PhotonBunch
        """
        if type(to_be_added) == PolarizedPhoton:
            self.photon_bunch.append(to_be_added)

        elif type(to_be_added) == list:
            self.photon_bunch.extend(to_be_added)

        elif type(to_be_added) == PhotonBunch:
            self.photon_bunch.extend(to_be_added.photon_bunch)

        else:
            raise TypeError("The photon(s) could not be added to the bunch!")

        self._set_dict()  # Update the array_dict attribute.

    def __len__(self):
        return len(self.photon_bunch)

    def __iter__(self):
        return iter(self.photon_bunch)

    def __getitem__(self, key):
        return self.photon_bunch[key]

    def _set_dict(self):
        """
        defines a dictionary containing information about the bunch.
        """
        self.array_dict = dict()
        energies = np.zeros(len(self))
        deviations = np.zeros(len(self))
        stokes = np.zeros([4, len(self)])
        polarization_degrees = np.zeros(len(self))
        i = 0

        for polarized_photon in self:
            energies[i] = polarized_photon.energy()  # Photon.energy()
            deviations[i] = polarized_photon.deviation()
            stokes[0, i] = polarized_photon.stokesVector().s0
            stokes[1, i] = polarized_photon.stokesVector().s1
            stokes[2, i] = polarized_photon.stokesVector().s2
            stokes[3, i] = polarized_photon.stokesVector().s3
            polarization_degrees[i] = polarized_photon.polarizationDegree()
            i += 1

        self.array_dict["number of photons"] = i
        self.array_dict["energies"] = energies
        self.array_dict["deviations"] = deviations
        self.array_dict["s0"] = stokes[0, :]
        self.array_dict["s1"] = stokes[1, :]
        self.array_dict["s2"] = stokes[2, :]
        self.array_dict["s3"] = stokes[3, :]
        self.array_dict["polarization degree"] = polarization_degrees

    def get_array(self, key):
        """
        :param key: 'deviations', 's0', 's1', 's2', 's3'.
        :return: numpy.ndarray
        """
        return self.array_dict[key]

    def is_monochromatic(self, places):
        """
        :param places: number of decimal places to be taken into account.
        :return: True if the bunch holds photons of the same energy.
        """
        first_energy = round(self.photon_bunch[0].energy(), places)

        # if the first element has the same energy as all others, then all others share the same energy value.
        for polarized_photon in self:
            if first_energy != round(polarized_photon.energy(), places):
                return False

        return True

    def is_unidirectional(self):
        """
        :return: True if the bunch holds photons going the same direction.
        """
        first_direction = self.photon_bunch[0].unitDirectionVector()  # Vector object.

        # if the first element goes the same direction as all others, then all others share the same direction.
        for polarized_photon in self:
            if first_direction != polarized_photon.unitDirectionVector():  # the precision is set to 7 decimal places.
                return False

        return True

    def to_string(self):
        """
        :return: string containing the parameters characterizing each photon in the bunch.
        """
        bunch_string = str()
        for polarized_photon in self:
            string_to_attach = str(polarized_photon.energy()) + " " + \
                               polarized_photon.unitDirectionVector().to_string() + " " + \
                               polarized_photon.stokesVector().to_string() + " " + \
                               str(polarized_photon.polarizationDegree()) + "\n"
            bunch_string += string_to_attach
        return bunch_string
