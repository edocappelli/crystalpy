"""
Represents a photon.
Except for energy all units are in SI. Energy is in eV.
"""
import scipy.constants.codata


class Photon(object):

    def __init__(self, energy_in_ev, direction_vector):
        """
        Constructor.
        :param energy_in_ev: Photon energy in eV.
        :param direction_vector: The direction of the photon as Vector.
        """
        self._energy_in_ev = float(energy_in_ev)
        self._unit_direction_vector = direction_vector.getNormalizedVector()

    def energy(self):
        """
        :return: Energy in eV.
        """
        return self._energy_in_ev

    def wavelength(self):
        """
        :return: The photon wavelength in meter.
        """
        codata = scipy.constants.codata.physical_constants
        speed_of_light = codata["speed of light in vacuum"][0]
        planck_constant = codata["Planck constant"][0]
        elementary_charge = codata["elementary charge"][0]
        E_in_Joule = self.energy() * elementary_charge

        # Wavelength in meter
        wavelength = (speed_of_light * planck_constant / E_in_Joule)

        return wavelength

    def wavenumber(self):
        """
        :return: Wavenumber in m^-1.
        """
        return (2.0 * scipy.constants.codata.pi) / self.wavelength()

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
