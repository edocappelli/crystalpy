"""
Unittest for CrystalPhasePlate class.
"""

import unittest

import numpy as np

from crystalpy.polarization.CrystalPhasePlate import CrystalPhasePlate
from crystalpy.util.StokesVector import StokesVector


class CrystalPhasePlateTest(unittest.TestCase):

    def test_create_matrix1(self):
        # phase plate 1

        element_list = [1, 1, 0, 0]
        incoming_stokes_vector = StokesVector(element_list)
        intensity_sigma = 0.5366453389170862
        phase_sigma = 1.9718541179322164  # radians
        intensity_pi = 0.9685465109320982
        phase_pi = -2.8903921744707217  # radians
        inclination_angle = 45.0 * np.pi / 180  # radians

        phase_plate = CrystalPhasePlate( #incoming_stokes_vector=incoming_stokes_vector,
                                        intensity_sigma=intensity_sigma, phase_sigma=phase_sigma,
                                        intensity_pi=intensity_pi, phase_pi=phase_pi,
                                        inclination_angle=inclination_angle)


        phase_plate_matrix = phase_plate.matrix


        self.assertIsInstance(phase_plate_matrix, np.ndarray)
        candidate = np.array([[0.7525959249245922, 0.0, -0.215950586007506, 0.0],
                             [-0.215950586007506, 0.0, 0.7525959249245922, 0.0],
                             [0.0, -0.10763540116609252, 0.0, -0.7128678636549213],
                             [0.0, -0.7128678636549213, 0.0, 0.10763540116609252]])
        # default precision: 6 decimal places.
        np.testing.assert_array_almost_equal(candidate, phase_plate_matrix)



    def test_calculate_stokes_vector(self):
        # phase plate 1
        # outgoing_stokes_vector = self.phase_plate1.calculate_stokes_vector()

        element_list = [1, 1, 0, 0]
        incoming_stokes_vector = StokesVector(element_list)
        intensity_sigma = 0.5366453389170862
        phase_sigma = 1.9718541179322164  # radians
        intensity_pi = 0.9685465109320982
        phase_pi = -2.8903921744707217  # radians
        inclination_angle = 45.0 * np.pi / 180  # radians

        phase_plate = CrystalPhasePlate( #incoming_stokes_vector=incoming_stokes_vector,
                                        intensity_sigma=intensity_sigma, phase_sigma=phase_sigma,
                                        intensity_pi=intensity_pi, phase_pi=phase_pi,
                                        inclination_angle=inclination_angle)



        outgoing_stokes_vector = phase_plate.calculate_stokes_vector(incoming_stokes_vector)

        self.assertIsInstance(outgoing_stokes_vector, StokesVector)
        candidate = np.array([0.7525959249245922,
                              -0.215950586007506,
                              -0.10763540116609252,
                              -0.7128678636549213])
        # default precision: 6 decimal places.
        np.testing.assert_array_almost_equal(candidate, outgoing_stokes_vector.components())
