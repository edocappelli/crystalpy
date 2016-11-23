"""
Unittest for CrystalPhasePlate class.
"""

import unittest

from orangecontrib.crystal.polarization.CrystalPhasePlate import CrystalPhasePlate
from orangecontrib.crystal.polarization.StokesVector import StokesVector

import numpy as np


def _generate_phase_plate1():

    element_list = [1, 1, 0, 0]
    incoming_stokes_vector = StokesVector(element_list)
    intensity_sigma = 0.5366453389170862
    phase_sigma = 1.9718541179322164  # radians
    intensity_pi = 0.9685465109320982
    phase_pi = -2.8903921744707217  # radians
    inclination_angle = 45.0  # degrees

    phase_plate = CrystalPhasePlate(incoming_stokes_vector=incoming_stokes_vector,
                                    intensity_sigma=intensity_sigma, phase_sigma=phase_sigma,
                                    intensity_pi=intensity_pi, phase_pi=phase_pi,
                                    inclination_angle=inclination_angle)
    return phase_plate


def _generate_phase_plate2():

    element_list = [-0.9716253979537108 * 1e-6, 0.6443153476193257 * 1e-6,
                    -0.45770218408141883 * 1e-6, -0.9302554106749403 * 1e-6]
    incoming_stokes_vector = StokesVector(element_list)
    intensity_sigma = 0.5366453389170862 * 1e-6
    phase_sigma = 1.9718541179322164 * 1e-6  # radians
    intensity_pi = 0.9685465109320982 * 1e-6
    phase_pi = -2.8903921744707217 * 1e-6  # radians
    inclination_angle = 90.0  # degrees

    phase_plate = CrystalPhasePlate(incoming_stokes_vector=incoming_stokes_vector,
                                    intensity_sigma=intensity_sigma, phase_sigma=phase_sigma,
                                    intensity_pi=intensity_pi, phase_pi=phase_pi,
                                    inclination_angle=inclination_angle)
    return phase_plate


class CrystalPhasePlateTest(unittest.TestCase):
    def setUp(self):
        self.phase_plate1 = _generate_phase_plate1()
        self.phase_plate2 = _generate_phase_plate2()

    def test_constructor(self):
        # phase plate 1
        self.assertListEqual(self.phase_plate1.incoming_stokes_vector.get_array(numpy=False), [1, 1, 0, 0])
        self.assertEqual(self.phase_plate1.intensity_sigma, 0.5366453389170862)
        self.assertEqual(self.phase_plate1.phase_sigma, 1.9718541179322164)
        self.assertEqual(self.phase_plate1.intensity_pi, 0.9685465109320982)
        self.assertEqual(self.phase_plate1.phase_pi, -2.8903921744707217)
        self.assertEqual(self.phase_plate1.inclination_angle, 45.0)
        # phase plate 2
        np.testing.assert_array_almost_equal(self.phase_plate2.incoming_stokes_vector.get_array(numpy=True),
                                             [-0.9716253979537108 * 1e-6,
                                              0.6443153476193257 * 1e-6,
                                              -0.45770218408141883 * 1e-6,
                                              -0.9302554106749403 * 1e-6])
        self.assertEqual(self.phase_plate2.intensity_sigma, 0.5366453389170862 * 1e-6)
        self.assertEqual(self.phase_plate2.phase_sigma, 1.9718541179322164 * 1e-6)
        self.assertEqual(self.phase_plate2.intensity_pi, 0.9685465109320982 * 1e-6)
        self.assertEqual(self.phase_plate2.phase_pi, -2.8903921744707217 * 1e-6)
        self.assertEqual(self.phase_plate2.inclination_angle, 90.0)

    def test_create_matrix(self):
        # phase plate 1
        phase_plate_matrix = self.phase_plate1._create_matrix()
        self.assertIsInstance(phase_plate_matrix, np.ndarray)
        candidate = np.array([[0.7525959249245922, 0.0, -0.215950586007506, 0.0],
                             [-0.215950586007506, 0.0, 0.7525959249245922, 0.0],
                             [0.0, -0.10763540116609252, 0.0, -0.7128678636549213],
                             [0.0, -0.7128678636549213, 0.0, 0.10763540116609252]])
        # default precision: 6 decimal places.
        np.testing.assert_array_almost_equal(candidate, phase_plate_matrix)

    def test_calculate_stokes_vector(self):
        # phase plate 1
        outgoing_stokes_vector = self.phase_plate1.calculate_stokes_vector()
        self.assertIsInstance(outgoing_stokes_vector, StokesVector)
        candidate = np.array([0.7525959249245922,
                              -0.215950586007506,
                              -0.10763540116609252,
                              -0.7128678636549213])
        # default precision: 6 decimal places.
        np.testing.assert_array_almost_equal(candidate, outgoing_stokes_vector.get_array(numpy=True))
