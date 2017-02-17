"""
Represents Mueller calculation results.
"""
import numpy

# TODO inheritate from DiffractionResults?
class MuellerResult(object):

    def __init__(self, diffraction_result):
        """
        Constructor.
        :param diffraction_result: result of the diffraction from our setup.
        """
        self.diffraction_result = diffraction_result
        self.diffraction_setup = diffraction_result.diffractionSetup()

        number_energies = len(self.energies())
        number_angles = len(self.angle_deviations())

        # Stokes parameters.
        self._s0 = numpy.zeros((number_energies,
                                number_angles))

        self._s1 = numpy.zeros((number_energies,
                                number_angles))

        self._s2 = numpy.zeros((number_energies,
                                number_angles))

        self._s3 = numpy.zeros((number_energies,
                                number_angles))

        # degree of circular polarization.
        self._circular_polarization_degree = numpy.zeros((number_energies,
                                                 number_angles))

    def energies(self):
        """
        Returns the energies used for these results.
        :return: Energies used for these results.
        """
        return self.diffraction_result.energies()

    def _energy_index(self, energy):
        """
        Returns the index of the entry in the energies list that is closest to the given energy.
        :param energy: Energy to find index for.
        :return: Energy index that corresponds to the energy.
        """
        energy_index = abs(self.energies()-energy).argmin()
        return energy_index

    def angle_deviations(self):
        """
        Returns the angle deviations used for these results.
        :return: Angle deviations used for these results.
        """
        return self.diffraction_result.angleDeviations()

    def _deviation_index(self, deviation):
        """
        Returns the index of the entry in the angle deviations list that is closest to the given deviation.
        :param deviation: Deviation to find index for.
        :return: Deviation index that corresponds to the deviation.
        """
        deviation_index = abs(self.angle_deviations()-deviation).argmin()
        return deviation_index

    def s0_by_energy(self, energy):
        """
        Returns the S0 Stokes parameter.
        :param energy: Energy to return S0 for.
        :return: S0.
        """
        energy_index = self._energy_index(energy)
        return self._s0[energy_index, :]

    def s1_by_energy(self, energy):
        """
        Returns the S1 Stokes parameter.
        :param energy: Energy to return S1 for.
        :return: S1.
        """
        energy_index = self._energy_index(energy)
        return self._s1[energy_index, :]

    def s2_by_energy(self, energy):
        """
        Returns the S2 Stokes parameter.
        :param energy: Energy to return S2 for.
        :return: S2.
        """
        energy_index = self._energy_index(energy)
        return self._s2[energy_index, :]

    def s3_by_energy(self, energy):
        """
        Returns the S3 Stokes parameter.
        :param energy: Energy to return S3 for.
        :return: S3.
        """
        energy_index = self._energy_index(energy)
        return self._s3[energy_index, :]

    def polarization_degree_by_energy(self, energy):
        """
        Returns the degree of circular polarization.
        :param energy: Energy to return the degree of circular polarization for.
        :return: degree of circular polarization.
        """
        energy_index = self._energy_index(energy)
        return self._circular_polarization_degree[energy_index, :]

    def s0_by_deviation(self, deviation):
        """
        Returns the S0 Stokes parameter.
        :param deviation: Deviation to return phase for.
        :return: S0.
        """
        deviation_index = self._deviation_index(deviation)
        return self._s0[deviation_index, :]

    def s1_by_deviation(self, deviation):
        """
        Returns the S1 Stokes parameter.
        :param deviation: Deviation to return phase for.
        :return: S1.
        """
        deviation_index = self._deviation_index(deviation)
        return self._s1[deviation_index, :]

    def s2_by_deviation(self, deviation):
        """
        Returns the S2 Stokes parameter.
        :param deviation: Deviation to return phase for.
        :return: S2.
        """
        deviation_index = self._deviation_index(deviation)
        return self._s2[deviation_index, :]

    def s3_by_deviation(self, deviation):
        """numpy.zeros((number_energies,
                                number_angles))
        Returns the S3 Stokes parameter.
        :param deviation: Deviation to return phase for.
        :return: S3.
        """
        deviation_index = self._deviation_index(deviation)
        return self._s3[deviation_index, :]

    def polarization_degree_by_deviation(self, deviation):
        """
        Returns the degree of circular polarization.
        :param deviation: Deviation to return the degree of circular polarization for.
        :return: degree of circular polarization.
        """
        deviation_index = self._deviation_index(deviation)
        return self._circular_polarization_degree[deviation_index, :]

    def add(self, energy, deviation, stokes_vector):
        """
        Adds a result for a given energy and deviation.
        """
        energy_index = self._energy_index(energy)
        deviation_index = self._deviation_index(deviation)

        self._s0[energy_index, deviation_index] = stokes_vector.s0
        self._s1[energy_index, deviation_index] = stokes_vector.s1
        self._s2[energy_index, deviation_index] = stokes_vector.s2
        self._s3[energy_index, deviation_index] = stokes_vector.s3
        self._circular_polarization_degree[energy_index, deviation_index] = stokes_vector.circularPolarizationDegree()
