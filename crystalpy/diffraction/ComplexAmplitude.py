"""
Represents a complex amplitude of a diffracted or transmissioned wave.
"""
import math


class ComplexAmplitude(object):

    def __init__(self, complex_amplitude):
        """
        Constructor.
        :param complex_amplitude: Complex amplitude of the wave.
        """
        self.setComplexAmplitude(complex_amplitude)

    def setComplexAmplitude(self, complex_amplitude):
        """
        Sets the complex amplitude.
        :param complex_amplitude: Complex amplitude of the wave.
        """
        self._complex_amplitude = complex_amplitude

    def rescale(self, scalar):
        """
        Rescales the complex amplitude.
        :param scalar: Scalar to rescale the complex amplitude with.
        """
        self._complex_amplitude = self._complex_amplitude * scalar

    def complexAmplitude(self):
        """
        Returns the complex amplitude.
        :return: Complex amplitude.
        """
        return self._complex_amplitude

    def intensity(self):
        """
        Return the intensity corresponding to the complex amplitude.
        :return: Intensity corresponding to the complex amplitude.
        """
        return abs(self._complex_amplitude) ** 2

    def phase(self):
        """
        Returns the phase of the complex amplitude.
        :return: Phase of the complex amplitude.
        """
        PP = self._complex_amplitude.real
        QQ = self._complex_amplitude.imag
        return math.atan2(QQ, PP)  # result between -pi and pi.

    def __truediv__(self, divisor):
        """
        Defines complex amplitude division.
        :param divisor: ComplexAmplitude dividing this instance.
        :return: Result of the division.
        """
        division = ComplexAmplitude(self.complexAmplitude() /
                                    divisor.complexAmplitude())
        return division
