"""
Represents a Mueller matrix.
See, e.g., https://en.wikipedia.org/wiki/Mueller_calculus
"""
import numpy

from crystalpy.util.StokesVector import StokesVector


class MuellerMatrix(object):

    def __init__(self, matrix=numpy.zeros((4,4)) ):
        """
        Constructor.
        :param matrix: Matrix as a numpy.ndarray object.
        """
        self.matrix = matrix

    @classmethod
    def initialize_as_general_linear_polarizer(cls,theta=0.0):
        mm = MuellerMatrix()
        mm.set_general_linear_polarizer(theta)
        return mm

    @classmethod
    def initialize_as_linear_polarizer_horizontal(cls):
        return cls.initialize_as_general_linear_polarizer(0.0)

    @classmethod
    def initialize_as_linear_polarizer_vertical(cls):
        return cls.initialize_as_general_linear_polarizer(numpy.pi/2)

    @classmethod
    def initialize_as_linear_polarizer_plus45(cls):
        return cls.initialize_as_general_linear_polarizer(numpy.pi/4)

    @classmethod
    def initialize_as_linear_polarizer_minus45(cls):
        return cls.initialize_as_general_linear_polarizer(-numpy.pi/4)

    @classmethod
    def initialize_as_general_linear_retarder(cls,theta=0.0, delta=0.0):
        mm = MuellerMatrix()
        mm.set_general_linear_retarder(theta,delta)
        return mm

    @classmethod
    def initialize_as_quarter_wave_plate_fast_vertical(cls):
        return cls.initialize_as_general_linear_retarder(numpy.pi/2,-numpy.pi/2)

    @classmethod
    def initialize_as_quarter_wave_plate_fast_horizontal(cls):
        return cls.initialize_as_general_linear_retarder(0.0,-numpy.pi/2)

    @classmethod
    def initialize_as_half_wave_plate(cls):
        return cls.initialize_as_general_linear_retarder(0.0,numpy.pi)

    @classmethod
    def initialize_as_ideal_mirror(cls):
        return cls.initialize_as_general_linear_retarder(0.0,numpy.pi)

    @classmethod
    def initialize_as_filter(cls,transmission=1.0):
        return cls.initialize_as_general_linear_retarder(0.0,0.0).matrix_by_scalar(transmission)


    def from_matrix_to_elements(self, return_numpy=True):
        """
        Returns a numpy.ndarray of the elements of a given matrix.
        If a list is needed one can use the numpy.array.tolist() method.
        :param numpy: if True returns numpy.ndarray, if False returns list.
        :return: [m00, m01, m02....mN0, mN1, mN2...]
        """
        matrix = numpy.asarray(self.matrix)
        result = matrix.flatten()

        if return_numpy:
            return result

        return list(result)

    def get_matrix(self):
        return self.matrix

    def matrix_by_scalar(self, scalar):
        """
        Multiplies the matrix by a scalar.
        :param scalar: the scalar factor.
        :return: new Mueller matrix.
        """
        new_mueller_matrix = self.matrix * scalar

        return MuellerMatrix(new_mueller_matrix)

    def matrix_by_vector(self, vector, return_numpy=True):
        """
        Multiplies the matrix by a vector.
        :param numpy: if True returns numpy.ndarray, if False returns list.
        :param vector: the vector factor.
        :return: matrix * vector (not a MuellerMatrix object).
        """
        matrix = numpy.asarray(self.matrix)
        result = numpy.dot(matrix, vector)

        if return_numpy:
            return result

        return list(result)

    def vector_by_matrix(self, vector, return_numpy=True):
        """
        Multiplies the matrix by a vector.
        :param numpy: if True returns numpy.ndarray, if False returns list.
        :param vector: the vector factor.
        :return: matrix * vector (not a MuellerMatrix object).
        """
        matrix = numpy.asarray(self.matrix)
        result = numpy.dot(vector, matrix)

        if return_numpy:
            return result

        return list(result)

    def mueller_times_mueller(self, matrix_2, mod=False):
        """
        Multiplies two Mueller matrices.
        :param matrix_2: Mueller matrix factor.
        :param mod: matrix multiplication is not commutative
                 -> mod controls which of the two matrices is the first factor.
        :return: Mueller matrix product.
        """
        matrix_1 = self.matrix

        if mod:
            product = numpy.dot(matrix_1, matrix_2)

        else:
            product = numpy.dot(matrix_2, matrix_1)

        return MuellerMatrix(product)

    def __eq__(self, candidate):
        """
        Determines whether two Mueller matrices are equal.
        :param candidate: Mueller matrix to compare to.
        :return: True if equal. False if not.
        """
        for i in range(4):
            for j in range(4):

                if self.matrix[i, j] != candidate.matrix[i, j]:
                    return False

        return True

    def __ne__(self, candidate):
        """
        Determines whether two Mueller matrices are not equal.
        :param candidate: Mueller matrix to compare to.
        :return: True if not equal. False if equal.
        """
        return not self == candidate

    def set_general_linear_polarizer(self,theta):

        """
        :param theta: angle of the polarizer in rad
        :return: Mueller matrix for a general liner polarizer (numpy.ndarray)
                See: https://en.wikipedia.org/wiki/Mueller_calculus
        """

        # First row.
        self.matrix[0, 0] = 0.5
        self.matrix[0, 1] = 0.5 * numpy.cos(2*theta)
        self.matrix[0, 2] = 0.5 * numpy.sin(2*theta)
        self.matrix[0, 3] = 0.0

        # Second row.
        self.matrix[1, 0] = 0.5  * numpy.cos(2*theta)
        self.matrix[1, 1] = 0.5  * (numpy.cos(2*theta))**2
        self.matrix[1, 2] = 0.5  * numpy.sin(2*theta) * numpy.cos(2*theta)
        self.matrix[1, 3] = 0.0

        # Third row.
        self.matrix[2, 0] = 0.5 * numpy.sin(2*theta)
        self.matrix[2, 1] = 0.5 * numpy.sin(2*theta) * numpy.cos(2*theta)
        self.matrix[2, 2] = 0.5 * (numpy.sin(2*theta))**2
        self.matrix[2, 3] = 0.0

        # Fourth row.
        self.matrix[3, 0] = 0.0
        self.matrix[3, 1] = 0.0
        self.matrix[3, 2] = 0.0
        self.matrix[3, 3] = 0.0


    def set_general_linear_retarder(self,theta,delta=0.0):

        """
        :param theta: angle of fast axis in rad
        :param delta: phase difference in rad between the fast and slow axis
        :return: Mueller matrix for a general liner polarizer (numpy.ndarray)
                See: https://en.wikipedia.org/wiki/Mueller_calculus
        """

        # First row.
        self.matrix[0, 0] = 1.0
        self.matrix[0, 1] = 0.0
        self.matrix[0, 2] = 0.0
        self.matrix[0, 3] = 0.0

        # Second row.    (numpy.cos(2*theta))**2  (numpy.sin(2*theta))**2 (numpy.cos(delta))**2 (numpy.sin(delta))**2
        self.matrix[1, 0] = 0.0
        self.matrix[1, 1] = (numpy.cos(2*theta))**2 + numpy.cos(delta) * (numpy.sin(2*theta))**2
        self.matrix[1, 2] = numpy.cos(2*theta) * numpy.sin(2*theta) - numpy.cos(2*theta) * numpy.cos(delta) * numpy.sin(2*theta)
        self.matrix[1, 3] = numpy.sin(2*theta) * numpy.sin(delta)

        # Third row.
        self.matrix[2, 0] = 0.0
        self.matrix[2, 1] =  numpy.cos(2*theta) * numpy.sin(2*theta) - numpy.cos(2*theta) * numpy.cos(delta) * numpy.sin(2*theta)
        self.matrix[2, 2] =  numpy.cos(delta) * (numpy.cos(2*theta))**2 + (numpy.sin(2*theta))**2
        self.matrix[2, 3] = -numpy.cos(2*theta) * numpy.sin(delta)

        # Fourth row.
        self.matrix[3, 0] = 0.0
        self.matrix[3, 1] = -numpy.sin(2*theta) * numpy.sin(delta)
        self.matrix[3, 2] =  numpy.cos(2*theta) * numpy.sin(delta)
        self.matrix[3, 3] =  numpy.cos(delta)

    def calculate_stokes_vector(self,incoming_stokes_vector):
        """
        Takes an incoming Stokes vector, multiplies it by a Mueller matrix
        and gives an outgoing Stokes vector as a result.
        :return: StokesVector object.
        """
        # incoming_stokes_vector = self.incoming_stokes_vector.get_array()  # Stokes vector.
        element_list = self.matrix_by_vector(incoming_stokes_vector.getList())
        outgoing_stokes_vector = StokesVector(element_list)

        return outgoing_stokes_vector
