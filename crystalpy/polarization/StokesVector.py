"""
Represents a Stokes vector.
Except for energy all units are in SI. Energy is in eV.
"""
from numpy import asarray


class StokesVector(object):

    def __init__(self, element_list):
        """
        Constructor.
        :param element_list: list containing the Stokes parameters S0,S1,S2,S3.
        """
        self.s0 = element_list[0]
        self.s1 = element_list[1]
        self.s2 = element_list[2]
        self.s3 = element_list[3]

    def get_array(self, numpy=True):
        """
        Generates a 1x4 array from the Stokes vector components.
        :param numpy: if True returns numpy.ndarray, if False returns list.
        :return: 1x4 array containing the Stokes parameters.
        """
        result = list()
        result.append(self.s0)
        result.append(self.s1)
        result.append(self.s2)
        result.append(self.s3)

        if numpy:
            return asarray(result)

        return result

    def polarization_degree(self):
        """
        Calculates the degree of circular polarization of the radiation
        described by the Stokes parameter.
        :return: degree of circular polarization
        """
        return self.s3 / self.s0

    def __eq__(self, candidate):
        """
        Determines whether two Stokes vectors are equal.
        :param candidate: Stokes vector to compare to.
        :return: True if equal. False if not.
        """
        if self.s0 != candidate.s0:
            return False

        if self.s1 != candidate.s1:
            return False

        if self.s2 != candidate.s2:
            return False

        if self.s3 != candidate.s3:
            return False

        return True

    def __ne__(self, candidate):
        """
        Determines whether two Stokes vectors are not equal.
        :param candidate: Stokes vector to compare to.
        :return: True if not equal. False if equal.
        """
        return not self == candidate
