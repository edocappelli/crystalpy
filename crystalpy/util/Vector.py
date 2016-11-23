"""
Represents a 3d vector.
"""
import numpy as np


class Vector(object):
    def __init__(self, x, y, z):
        """
        Constructor.
        :param x: x component.
        :param y: y component.
        :param z: z component.
        """
        self.setComponents(x, y, z)

    @staticmethod
    def fromComponents(components):
        """
        Creates a vector from a list/array of at least three elements.
        :param components: x,y,z components of the vector.
        :return: Vector having these x,y,z components.
        """
        return Vector(components[0],
                      components[1],
                      components[2])

    def setComponents(self, x, y, z):
        """
        Sets vector components.
        :param x: x component.
        :param y: y component.
        :param z: z component.
        """
        self._components = np.asarray([x, y, z])

    def components(self):
        """
        Returns the components of this vector a 1d three element array.
        :return:
        """
        return self._components

    def __eq__(self, candidate):
        """
        Determines if two vectors are equal.
        :param candidate: Vector to compare to.
        :return: True if both vectors are equal. Otherwise False.
        """
        return np.linalg.norm(self.components()
                              -
                              candidate.components()) < 1.e-7

    def __ne__(self, candidate):
        """
        Determines if two vectors are not equal.
        :param candidate: Vector to compare to.
        :return: True if both vectors are not equal. Otherwise False.
        """
        return not (self == candidate)

    def addVector(self, summand):
        """
        Adds two vectors.
        :param summand: The vector to add to this instance.
        :return: The sum as a vector.
        """
        components = self.components() + summand.components()
        return Vector.fromComponents(components)

    def scalarMultiplication(self, factor):
        """
        Scalar multiplies this vector.
        :param factor: The scalar to multiply with.
        :return: Scalar multiplied vector.
        """
        components = self.components() * factor
        return Vector.fromComponents(components)

    def subtractVector(self, subtrahend):
        """
        Subtract a vector from this instance.
        :param subtrahend: Vector to subtract.
        :return: The difference of the two vectors.
        """
        result = self.addVector(subtrahend.scalarMultiplication(-1.0))
        return result

    def scalarProduct(self, factor):
        """
        Calculates the scalar product of this vector with the given vector.
        :param factor: The vector to calculate the scalar product with.
        :return: Scalar product of the two vectors.
        """
        scalar_product = np.dot(self.components(), factor.components())
        return scalar_product

    def crossProduct(self, factor):
        """
        Calculates the cross product of two vectors.
        :param factor: The vector to form the cross product with.
        :return: Cross product of the two vectors.
        """
        components = np.cross(self.components(), factor.components())
        return Vector.fromComponents(components)

    def norm(self):
        """
        Returns the standard norm of this norm.
        :return: Norm of this vector,
        """
        norm = self.scalarProduct(self) ** 0.5
        return norm

    def getNormalizedVector(self):
        """
        Returns a normalized vector of this vector.
        :return: Normalized vector of this vector.
        """
        return self.scalarMultiplication(self.norm() ** -1.0)

    def rotateAroundAxis(self, rotation_axis, angle):
        """
        Rotates the vector around an axis.
        :param rotation_axis: Vector specifying the rotation axis.
        :param angle: Rotation angle.
        :return: Rotated vector.
        """
        # For the mathematics look for: Rodrigues rotation formula.
        # http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        unit_rotation_axis = rotation_axis.getNormalizedVector()

        rotated_vector = self.scalarMultiplication(np.cos(angle))

        tmp_vector = unit_rotation_axis.crossProduct(self)
        tmp_vector = tmp_vector.scalarMultiplication(np.sin(angle))
        rotated_vector = rotated_vector.addVector(tmp_vector)

        scalar_factor = self.scalarProduct(unit_rotation_axis) * (1.0 - np.cos(angle))
        tmp_vector = unit_rotation_axis.scalarMultiplication(scalar_factor)
        rotated_vector = rotated_vector.addVector(tmp_vector)

        return rotated_vector

    def parallelTo(self, vector):
        """
        Returns the parallel projection of this vector along the given vector.
        :param vector: Vector defining the parallel direction.
        :return: Parallel projection along the vector.
        """
        unit_direction = vector.getNormalizedVector()
        projection_in_direction = self.scalarProduct(unit_direction)
        parallel_projection = unit_direction.scalarMultiplication(projection_in_direction)

        return parallel_projection

    def perpendicularTo(self, vector):
        """
        Returns the projection perpendicular to the given vector.
        :param vector: Vector that defines the direction.
        :return: Projection perpendicular to the given vector.
        """
        perpendicular = self.subtractVector(self.parallelTo(vector))
        return perpendicular

    def getOnePerpendicularVector(self):
        """
        Returns one arbitrary vector perpendicular to this vector.
        :return: One arbitrary vector perpendicular to this vector.
        """
        vector_y = Vector(0, 1, 0)
        vector_z = Vector(0, 0, 1)

        if self.getNormalizedVector() == vector_z:
            return vector_y

        vector_perpendicular = vector_z.perpendicularTo(self)
        vector_perpendicular = vector_perpendicular.getNormalizedVector()

        return vector_perpendicular

    def angle(self, factor):
        """
        Return the angle between this vector and the given vector.
        :param factor: Vector to determine the angle with.
        :return: Angle between this vector and the given vector.
        """
        n1 = self.getNormalizedVector()
        n2 = factor.getNormalizedVector()

        # Determine angle between the two vectors.
        cos_angle = n1.scalarProduct(n2)
        angle = np.arccos(cos_angle)
        # Edoardo: numpy.arccos() always returns an angle in radians in [0, pi].

        # Mark's version:
        # By convention always return the smaller angle.
        # while angle > 2.0 * np.pi:
            # angle -= 2.0 * np.pi

        # if angle > np.pi:
            # angle = 2.0 * np.pi - angle

        return angle

    def getVectorWithAngle(self, angle):
        """
        Returns one arbitrary vector with the given angle.
        :param angle: The requested angle.
        :return:Vector with given angle to this vector.
        """
        vector_perpendicular = self.getOnePerpendicularVector()
        vector_with_angle = self.rotateAroundAxis(vector_perpendicular,
                                                  angle)

        return vector_with_angle

    def printComponents(self):
        """
        Prints the components of this vector.
        """
        print("Vector: x", self.components()[0])
        print("Vector: y", self.components()[1])
        print("Vector: z", self.components()[2])
