from crystalpy.polarization.MuellerMatrix import MuellerMatrix
import numpy


# TODO rename to PerfectCrystalMuellerMatrix
class CrystalPhasePlate(MuellerMatrix):

    def __init__(self, # incoming_stokes_vector,
                 intensity_sigma, phase_sigma,
                 intensity_pi, phase_pi,
                 inclination_angle=0.0):
        """
        Constructor.
        """

        # self.incoming_stokes_vector = incoming_stokes_vector  # StokesVector object.

        phase_plate_matrix = self._create_matrix(intensity_sigma, phase_sigma,
                 intensity_pi, phase_pi,
                 inclination_angle)
        super(CrystalPhasePlate, self).__init__(phase_plate_matrix)

    def _create_matrix(self,
                        intensity_sigma, phase_sigma,
                        intensity_pi, phase_pi,
                        inclination_angle):
        """
        TODO: put article with the notation
        :return: Mueller matrix for a phase plate (numpy.ndarray).
        """
        alpha = inclination_angle  # radians.

        # Create the Mueller matrix for a phase plate as a numpy array.
        phase_plate_matrix = numpy.zeros([4, 4])

        # First row.
        phase_plate_matrix[0, 0] = 0.5 * (intensity_sigma + intensity_pi)
        phase_plate_matrix[0, 1] = 0.5 * (intensity_sigma - intensity_pi) * numpy.cos(2 * alpha)
        phase_plate_matrix[0, 2] = 0.5 * (intensity_sigma - intensity_pi) * numpy.sin(2 * alpha)
        phase_plate_matrix[0, 3] = 0.0

        # Second row.
        phase_plate_matrix[1, 0] = 0.5 * (intensity_sigma - intensity_pi)
        phase_plate_matrix[1, 1] = 0.5 * (intensity_sigma + intensity_pi) * numpy.cos(2 * alpha)
        phase_plate_matrix[1, 2] = 0.5 * (intensity_sigma + intensity_pi) * numpy.sin(2 * alpha)
        phase_plate_matrix[1, 3] = 0.0

        scalar = numpy.sqrt(intensity_sigma) * numpy.sqrt(intensity_pi)
        delta_phase = phase_pi - phase_sigma

        # Third row.
        phase_plate_matrix[2, 0] = 0.0
        phase_plate_matrix[2, 1] = - scalar * numpy.cos(delta_phase) * numpy.sin(2 * alpha)
        phase_plate_matrix[2, 2] = scalar * numpy.cos(delta_phase) * numpy.cos(2 * alpha)
        phase_plate_matrix[2, 3] = - scalar * numpy.sin(delta_phase)

        # Fourth row.
        phase_plate_matrix[3, 0] = 0.0
        phase_plate_matrix[3, 1] = - scalar * numpy.sin(delta_phase) * numpy.sin(2 * alpha)
        phase_plate_matrix[3, 2] = scalar * numpy.sin(delta_phase) * numpy.cos(2 * alpha)
        phase_plate_matrix[3, 3] = scalar * numpy.cos(delta_phase)

        return phase_plate_matrix

