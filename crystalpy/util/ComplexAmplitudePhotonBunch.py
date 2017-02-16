"""
This object contains a list of PolarizedPhoton objects, characterized by energy, direction vector and Stokes vector.
This object is used as input to and output from the passive crystal widget.
"""
import numpy as np

from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitude
from crystalpy.util.PhotonBunch import PhotonBunch


class ComplexAmplitudePhotonBunch(PhotonBunch):
    """
    is a collection of ComplexAmplitudePhoton objects, making up the photon beam.
    """
    def __init__(self, complex_amplitude_photons=None):
        """
        :param polarized_photons: bunch of PolarizedPhoton objects.
        :type polarized_photons: list(PolarizedPhoton, PolarizedPhoton, ...)
        """
        if complex_amplitude_photons == None:
            self.polarized_photon_bunch = []
        else:
            self.polarized_photon_bunch = complex_amplitude_photons


    def toDictionary(self):
        """
        defines a dictionary containing information about the bunch.
        """
        array_dict = PhotonBunch.toDictionary(self)

        intensityS = np.zeros(len(self))
        intensityP = np.zeros_like(intensityS)
        phaseS     = np.zeros_like(intensityS)
        phaseP     = np.zeros_like(intensityS)


        for i,polarized_photon in enumerate(self):
            intensityS[i] = polarized_photon.getIntensityS()
            intensityP[i] = polarized_photon.getIntensityP()
            phaseS    [i] = polarized_photon.getPhaseS()
            phaseP    [i] = polarized_photon.getPhaseP()

        array_dict["intensityS"] = intensityS
        array_dict["intensityP"] = intensityP
        array_dict["phaseS"] = phaseS
        array_dict["phaseP"] = phaseP


        return array_dict


    def toString(self):
        """
        :return: string containing the parameters characterizing each photon in the bunch.
        """
        bunch_string = str()

        for i in range(self.getNumberOfPhotons()):
            photon = self.getPhotonIndex(i)
            string_to_attach = str(photon.energy()) + " " + \
                               photon.unitDirectionVector().toString() + "\n"
            bunch_string += string_to_attach
        return bunch_string