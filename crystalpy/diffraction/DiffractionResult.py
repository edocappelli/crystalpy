"""
Represents diffraction results.
"""
import numpy


class DiffractionResult(object):

    INDEX_POLARIZATION_S = 0
    INDEX_POLARIZATION_P = 1
    INDEX_DIFFERENCE_SP = 2

    def __init__(self, diffraction_setup, bragg_angle):
        """
        Constructor.
        :param diffraction_setup: Setup used for these results.
        :param bragg_angle: Bragg angle of the setup.
        """
        self._diffraction_setup = diffraction_setup.clone()
        self._bragg_angle = bragg_angle

        number_energies = len(self.energies())
        number_angles = len(self.angleDeviations())
        number_polarizations = 3

        self._intensities = numpy.zeros((number_energies,
                                         number_angles,
                                         number_polarizations))

        self._phases = numpy.zeros((number_energies,
                                    number_angles,
                                    number_polarizations))

    def diffractionSetup(self):
        """
        Returns the diffraction setup used for the calculation of these results.
        :return: Diffraction setup used for the calculation of these results.
        """
        return self._diffraction_setup

    def braggAngle(self):
        """
        Returns Bragg angle used for these results.
        :return: Bragg angle used for these results.
        """
        return self._bragg_angle

    def energies(self):
        """
        Returns the energies used for these results.
        :return: Energies used for these results.
        """
        return self._diffraction_setup.energies()

    def _energyIndexByEnergy(self, energy):
        """
        Returns the index of the entry in the energies list that is closest to the given energy.
        :param energy: Energy to find index for.
        :return: Energy index that corresponds to the energy.
        """
        energy_index = abs(self.energies()-energy).argmin()
        return energy_index

    def angleDeviations(self):
        """
        Returns the angle deviations used for these results.
        :return: Angle deviations used for these results.
        """
        return self._diffraction_setup.angleDeviationGrid()

    def _deviationIndexByDeviation(self, deviation):
        """
        Returns the index of the entry in the angle deviations list that is closest to the given deviation.
        :param deviation: Deviation to find index for.
        :return: Deviation index that corresponds to the deviation.
        """
        deviation_index = abs(self.angleDeviations()-deviation).argmin()
        return deviation_index

    def angles(self):
        """
        Returns the angles used for calculation of these results.
        :return: The angles used for the calculation of these results.
        """
        return [self.braggAngle() + dev for dev in self.angleDeviations()]

    def sIntensityByEnergy(self, energy):
        """
        Returns the intensity of the S polarization.
        :param energy: Energy to return intensity for.
        :return: Intensity of the S polarization.
        """
        energy_index = self._energyIndexByEnergy(energy)
        return self._intensities[energy_index, :, self.INDEX_POLARIZATION_S]

    def sPhaseByEnergy(self, energy):
        """
        Returns the phase of the S polarization.
        :param energy: Energy to return phase for.
        :return: Phase of the S polarization.
        """
        energy_index = self._energyIndexByEnergy(energy)
        return self._phases[energy_index, :, self.INDEX_POLARIZATION_S]

    def pIntensityByEnergy(self, energy):
        """
        Returns the intensity of the P polarization.
        :param energy: Energy to return intensity for.
        :return: Intensity of the P polarization.
        """
        energy_index = self._energyIndexByEnergy(energy)
        return self._intensities[energy_index, :, self.INDEX_POLARIZATION_P]

    def pPhaseByEnergy(self, energy):
        """
        Returns the phase of the P polarization.
        :param energy: Energy to return phase for.
        :return: Phase of the P polarization.
        """
        energy_index = self._energyIndexByEnergy(energy)
        return self._phases[energy_index, :, self.INDEX_POLARIZATION_P]

    def differenceIntensityByEnergy(self, energy):
        """
        Returns the intensity of the difference between S and P polarizations.
        :param energy: Energy to return intensity for.
        :return: Intensity of the  difference between the S and P polarization.
        """
        energy_index = self._energyIndexByEnergy(energy)
        return self._intensities[energy_index, :, self.INDEX_DIFFERENCE_SP]

    def differencePhaseByEnergy(self, energy):
        """
        Returns the phase of the difference between S and P polarizations.
        :param energy: Energy to return phase for.
        :return: Phase of the difference between S and P polarization.
        """
        energy_index = self._energyIndexByEnergy(energy)
        return self._phases[energy_index, :, self.INDEX_DIFFERENCE_SP]

    def sIntensityByDeviation(self, deviation):
        """
        Returns the intensity of the S polarization.
        :param deviation: Deviation to return intensity for.
        :return: Intensity of the S polarization.
        """
        deviation_index = self._deviationIndexByDeviation(deviation)
        return self._intensities[:, deviation_index, self.INDEX_POLARIZATION_S]

    def sPhaseByDeviation(self, deviation):
        """
        Returns the phase of the S polarization.
        :param deviation: Deviation to return phase for.
        :return: Phase of the S polarization.
        """
        deviation_index = self._deviationIndexByDeviation(deviation)
        return self._phases[:, deviation_index, self.INDEX_POLARIZATION_S]

    def pIntensityByDeviation(self, deviation):
        """
        Returns the intensity of the P polarization.
        :param deviation: Deviation to return intensity for.
        :return: Intensity of the P polarization.
        """
        deviation_index = self._deviationIndexByDeviation(deviation)
        return self._intensities[:, deviation_index, self.INDEX_POLARIZATION_P]

    def pPhaseByDeviation(self, deviation):
        """
        Returns the phase of the P polarization.
        :param deviation: Deviation to return phase for.
        :return: Phase of the P polarization.
        """
        deviation_index = self._deviationIndexByDeviation(deviation)
        return self._phases[:, deviation_index, self.INDEX_POLARIZATION_P]

    def differenceIntensityByDeviation(self, deviation):
        """
        Returns the intensity of the difference between S and P polarizations.
        :param deviation: Deviation to return intensity for.
        :return: Intensity of the  difference between the S and P polarization.
        """
        deviation_index = self._deviationIndexByDeviation(deviation)
        return self._intensities[:, deviation_index, self.INDEX_DIFFERENCE_SP]

    def differencePhaseByDeviation(self, deviation):
        """
        Returns the phase of the difference between S and P polarizations.
        :param deviation: Deviation to return phase for.
        :return: Phase of the difference between S and P polarization.
        """
        deviation_index = self._deviationIndexByDeviation(deviation)
        return self._phases[:, deviation_index, self.INDEX_DIFFERENCE_SP]

    def add(self, energy, deviation, s_complex_amplitude, p_complex_amplitude, difference_complex_amplitude):
        """
        Adds a result for a given energy and deviation.
        """
        energy_index = self._energyIndexByEnergy(energy)
        deviation_index = self._deviationIndexByDeviation(deviation)

        self._intensities[energy_index, deviation_index, self.INDEX_POLARIZATION_S] = s_complex_amplitude.intensity()
        self._intensities[energy_index, deviation_index, self.INDEX_POLARIZATION_P] = p_complex_amplitude.intensity()
        self._intensities[energy_index, deviation_index, self.INDEX_DIFFERENCE_SP] = difference_complex_amplitude.intensity()

        self._phases[energy_index, deviation_index, self.INDEX_POLARIZATION_S] = s_complex_amplitude.phase()
        self._phases[energy_index, deviation_index, self.INDEX_POLARIZATION_P] = p_complex_amplitude.phase()
        self._phases[energy_index, deviation_index, self.INDEX_DIFFERENCE_SP] = difference_complex_amplitude.phase()

    # def _debugPlot(self):
        # """
        # Debug plot intensities.
        # """
        # x = [i * 1e+6 for i in self.angleDeviations()]
        # plot(x, self.sIntensityByEnergy(), label="S polarization")
        # plot(x, self.pIntensityByEnergy(), label="P polarization")
        # plot(x, self.sPhase(), label="S polarization")
        # plot(x, self.pPhase(), label="P polarization")
        # show()
        # return

    # def _debugPrint(self):
        # """
        # Debug print of S and P intensities and phases.
        # """
        # S polarization.
        # print("s_intensity="+str(self.sIntensityByEnergy()).replace("array(","").replace(") * dimensionless",""))
        # print("s_phase="+str(self.sPhaseByEnergy()))

        # P polarization.
        # print("p_intensity="+str(self.pIntensityByEnergy()).replace("array(","").replace(") * dimensionless",""))
        # print("p_phase="+str(self.pPhaseByEnergy()))
