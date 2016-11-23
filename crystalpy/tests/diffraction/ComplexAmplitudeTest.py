"""
Unittest for ComplexAmplitude class.
"""

import unittest

import numpy as np

from orangecontrib.crystal.diffraction.ComplexAmplitude import ComplexAmplitude


class ComplexAmplitudeTest(unittest.TestCase):
    def testConstructor(self):
        number = 1 + 2j
        complex_amplitude = ComplexAmplitude(number)
        
        self.assertAlmostEqual(complex_amplitude._complex_amplitude,
                               number)
        
    def testSetComplexAmplitude(self):
        complex_amplitude = ComplexAmplitude(0)
        number = 1 + 2j
        complex_amplitude.setComplexAmplitude(number)
        
        self.assertAlmostEqual(complex_amplitude._complex_amplitude,
                               number)

    def testIntensity(self):
        number = 1 + 2j
        complex_amplitude = ComplexAmplitude(number)

        self.assertAlmostEqual(complex_amplitude.intensity(),
                               5)

    def testRescale(self):
        number = 1 + 2j
        complex_amplitude = ComplexAmplitude(number)

        complex_amplitude.rescale(2.0)

        self.assertAlmostEqual(complex_amplitude._complex_amplitude,
                               2 + 4j)

    def testPhase(self):
        complex_amplitude = ComplexAmplitude(1 + 1j)

        self.assertAlmostEqual(complex_amplitude.phase(),
                               np.pi / 4.0)

    def testDivision(self):
        number_one = 1 + 2j
        complex_amplitude_one = ComplexAmplitude(number_one)

        number_two = 1 + 1j
        complex_amplitude_two = ComplexAmplitude(number_two)

        result = complex_amplitude_one / complex_amplitude_two

        self.assertAlmostEqual(result.intensity(), 2.5)
        self.assertAlmostEqual(result.phase(), 0.32175055)
