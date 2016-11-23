"""
Represents a Mueller matrix.
Except for energy all units are in SI. Energy is in eV.
"""
import numpy as np


class MuellerMatrix(object):

    def __init__(self, matrix):
        """
        Constructor.
        :param matrix: Matrix as a numpy.ndarray object.
        """
        self.matrix = matrix

    def from_matrix_to_elements(self, numpy=True):
        """
        Returns a numpy.ndarray of the elements of a given matrix.
        If a list is needed one can use the numpy.array.tolist() method.
        :param numpy: if True returns numpy.ndarray, if False returns list.
        :return: [m00, m01, m02....mN0, mN1, mN2...]
        """
        matrix = np.asarray(self.matrix)
        result = matrix.flatten()

        if numpy:
            return result

        return list(result)

    def matrix_by_scalar(self, scalar):
        """
        Multiplies the matrix by a scalar.
        :param scalar: the scalar factor.
        :return: new Mueller matrix.
        """
        new_mueller_matrix = self.matrix * scalar

        return MuellerMatrix(new_mueller_matrix)

    def matrix_by_vector(self, vector, numpy=True):
        """
        Multiplies the matrix by a vector.
        :param numpy: if True returns numpy.ndarray, if False returns list.
        :param vector: the vector factor.
        :return: matrix * vector (not a MuellerMatrix object).
        """
        matrix = np.asarray(self.matrix)
        result = np.dot(matrix, vector)

        if numpy:
            return result

        return list(result)

    def vector_by_matrix(self, vector, numpy=True):
        """
        Multiplies the matrix by a vector.
        :param numpy: if True returns numpy.ndarray, if False returns list.
        :param vector: the vector factor.
        :return: matrix * vector (not a MuellerMatrix object).
        """
        matrix = np.asarray(self.matrix)
        result = np.dot(vector, matrix)

        if numpy:
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
            product = np.dot(matrix_1, matrix_2)

        else:
            product = np.dot(matrix_2, matrix_1)

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
