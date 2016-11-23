from orangecontrib.crystal.polarization.StokesVector import StokesVector
from orangecontrib.crystal.polarization.MuellerMatrix import MuellerMatrix
import numpy as np


class CrystalPhasePlate(MuellerMatrix):

    def __init__(self, incoming_stokes_vector,
                 intensity_sigma, phase_sigma,
                 intensity_pi, phase_pi,
                 inclination_angle=0.0):
        """
        Constructor.
        :param intensity_sigma:
        :param phase_sigma:
        :param intensity_pi:
        :param phase_pi:
        :param inclination_angle:
        """
        self.intensity_sigma = intensity_sigma
        self.phase_sigma = phase_sigma
        self.intensity_pi = intensity_pi
        self.phase_pi = phase_pi
        self.inclination_angle = inclination_angle  # degrees.
        self.incoming_stokes_vector = incoming_stokes_vector  # StokesVector object.

        phase_plate_matrix = self._create_matrix()
        super(CrystalPhasePlate, self).__init__(phase_plate_matrix)

    def _create_matrix(self):
        """
        TODO: put article with the notation
        :return: Mueller matrix for a phase plate (numpy.ndarray).
        """
        alpha = self.inclination_angle * np.pi / 180  # degrees -> radians.

        # Create the Mueller matrix for a phase plate as a numpy array.
        phase_plate_matrix = np.zeros([4, 4])

        # First row.
        phase_plate_matrix[0, 0] = 0.5 * (self.intensity_sigma + self.intensity_pi)
        phase_plate_matrix[0, 1] = 0.5 * (self.intensity_sigma - self.intensity_pi) * np.cos(2 * alpha)
        phase_plate_matrix[0, 2] = 0.5 * (self.intensity_sigma - self.intensity_pi) * np.sin(2 * alpha)
        phase_plate_matrix[0, 3] = 0.0

        # Second row.
        phase_plate_matrix[1, 0] = 0.5 * (self.intensity_sigma - self.intensity_pi)
        phase_plate_matrix[1, 1] = 0.5 * (self.intensity_sigma + self.intensity_pi) * np.cos(2 * alpha)
        phase_plate_matrix[1, 2] = 0.5 * (self.intensity_sigma + self.intensity_pi) * np.sin(2 * alpha)
        phase_plate_matrix[1, 3] = 0.0

        scalar = np.sqrt(self.intensity_sigma) * np.sqrt(self.intensity_pi)
        delta_phase = self.phase_pi - self.phase_sigma

        # Third row.
        phase_plate_matrix[2, 0] = 0.0
        phase_plate_matrix[2, 1] = - scalar * np.cos(delta_phase) * np.sin(2 * alpha)
        phase_plate_matrix[2, 2] = scalar * np.cos(delta_phase) * np.cos(2 * alpha)
        phase_plate_matrix[2, 3] = - scalar * np.sin(delta_phase)

        # Fourth row.
        phase_plate_matrix[3, 0] = 0.0
        phase_plate_matrix[3, 1] = - scalar * np.sin(delta_phase) * np.sin(2 * alpha)
        phase_plate_matrix[3, 2] = scalar * np.sin(delta_phase) * np.cos(2 * alpha)
        phase_plate_matrix[3, 3] = scalar * np.cos(delta_phase)

        return phase_plate_matrix

    def calculate_stokes_vector(self):
        """
        Takes an incoming Stokes vector, multiplies it by a Mueller matrix
        and gives an outgoing Stokes vector as a result.
        :return: StokesVector object.
        """
        incoming_stokes_vector = self.incoming_stokes_vector.get_array()  # Stokes vector.
        element_list = self.matrix_by_vector(incoming_stokes_vector, numpy=False)
        outgoing_stokes_vector = StokesVector(element_list)

        return outgoing_stokes_vector
