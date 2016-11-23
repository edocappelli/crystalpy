"""
Unittest for Vector class.
"""
import unittest
from numpy import pi, arange, array

from orangecontrib.crystal.util.Vector import Vector


class VectorTest(unittest.TestCase):
    def testConstructor(self):
        vector = Vector(1, 2, 3)
        self.assertIsInstance(vector, Vector)

        self.assertAlmostEqual(vector.components()[0],
                               1)
        self.assertAlmostEqual(vector.components()[1],
                               2)
        self.assertAlmostEqual(vector.components()[2],
                               3)

    def testFromComponents(self):
        vector = Vector.fromComponents(array([11, -2, 23]))

        self.assertAlmostEqual(vector.components()[0],
                               11)
        self.assertAlmostEqual(vector.components()[1],
                               -2)
        self.assertAlmostEqual(vector.components()[2],
                               23)

    def testSetComponents(self):
        vector = Vector(1, 2, 3)

        vector.setComponents(-3.0, -4.0, 0.0)

        self.assertAlmostEqual(vector.components()[0],
                               -3.0)
        self.assertAlmostEqual(vector.components()[1],
                               -4.0)
        self.assertAlmostEqual(vector.components()[2],
                               0.0)

    def testComponents(self):
        vector = Vector(1, 2, 3)

        vector.setComponents(-3.0, -4.0, 0.0)

        self.assertAlmostEqual(vector.components()[0],
                               -3.0)
        self.assertAlmostEqual(vector.components()[1],
                               -4.0)
        self.assertAlmostEqual(vector.components()[2],
                               0.0)

    def testOperatorEqual(self):
        vector_1 = Vector(1, 2, 3)
        vector_2 = Vector(1.0, 2.0, 3.0)
        vector_3 = Vector(-1, 0, 0)

        self.assertTrue(vector_1 == vector_2)
        self.assertFalse(vector_1 == vector_3)

    def testOperatorNotEqual(self):
        vector_1 = Vector(1, 2, 3)
        vector_2 = Vector(1.0, 2.0, 3.0)
        vector_3 = Vector(-1, 0, 0)

        self.assertFalse(vector_1 != vector_2)
        self.assertTrue(vector_1 != vector_3)

    def testAddVector(self):
        vector_1 = Vector(1, 2, 3)
        vector_2 = Vector(-1, 2, 1)
        vector_sum = Vector(0, 4, 4)

        result = vector_1.addVector(vector_2)

        self.assertTrue(result == vector_sum)

    def testScalarMultiplication(self):
        vector_1 = Vector(1, 2, 3)
        vector_2 = Vector(3, 6, 9)
        result = Vector(-1, -2, -3)

        vector_1 = vector_1.scalarMultiplication(-1.0)

        self.assertTrue(result == vector_1)
        vector_1 = vector_1.scalarMultiplication(-3.0)
        self.assertTrue(vector_1 == vector_2)

    def testSubtractVector(self):
        vector_1 = Vector(1, 2, 3)
        vector_2 = Vector(-1, 2, 1)
        vector_diff = Vector(2, 0, 2)

        result = vector_1.subtractVector(vector_2)

        self.assertTrue(result == vector_diff)

    def testScalarProduct(self):
        vector_1 = Vector(1, 2, 3)
        vector_2 = Vector(-1, 2, 1)
        vector_3 = Vector(-1, 2, -1)

        result = vector_1.scalarProduct(vector_2)
        self.assertAlmostEqual(result,
                               6)

        result = vector_1.scalarProduct(vector_3)
        self.assertAlmostEqual(result,
                               0)

    def testCrossProduct(self):
        vector_x = Vector(1, 0, 0)
        vector_y = Vector(0, 1, 0)
        vector_z = Vector(0, 0, 1)
        vector_minus_z = Vector(0, 0, -1)
        vector_null = Vector(0, 0, 0)

        result = vector_x.crossProduct(vector_y)
        self.assertTrue(result == vector_z)

        result = vector_z.crossProduct(vector_x)
        self.assertTrue(result == vector_y)

        result = vector_y.crossProduct(vector_x)
        self.assertTrue(result == vector_minus_z)

        result = vector_x.crossProduct(vector_x)
        self.assertTrue(result == vector_null)

    def testNorm(self,):
        vector_6 = Vector(2, 4, 4)
        vector_2 = Vector(2, 0, 0)

        self.assertAlmostEqual(vector_6.norm(), 6.0)
        self.assertAlmostEqual(vector_2.norm(), 2.0)

    def testGetNormalizedVector(self):
        vector_6 = Vector(2, 4, 4)

        normalized_vector = vector_6.getNormalizedVector()
        self.assertAlmostEqual(normalized_vector.norm(), 1.0)

        self.assertTrue(vector_6 == normalized_vector.scalarMultiplication(6))

    def testRotateAroundAxis(self):
        vector_x = Vector(1, 0, 0)
        vector_y = Vector(0, 1, 0)
        vector_z = Vector(0, 0, 1)
        vector_minus_x = Vector(-1, 0, 0)
        vector_minus_y = Vector(0, -1, 0)
        vector_minus_z = Vector(0, 0, -1)


        result = vector_x.rotateAroundAxis(vector_z, pi)
        self.assertTrue(result == vector_minus_x)

        result = vector_x.rotateAroundAxis(vector_z, pi / 2.0)
        self.assertTrue(result == vector_y)

        result = vector_x.addVector(vector_y).rotateAroundAxis(vector_z, pi)
        self.assertTrue(result == vector_minus_x.addVector(vector_minus_y))

        result = vector_z.rotateAroundAxis(vector_z, pi)
        self.assertTrue(result == vector_z)

        result = vector_z.rotateAroundAxis(vector_x, pi / 2.0)
        self.assertTrue(result == vector_minus_y)

    def testParallelTo(self):
        vector = Vector(1, 1, 3)
        vector_z = Vector(0, 0, 1)

        result = vector.parallelTo(vector_z)

        self.assertTrue(result == Vector(0, 0, 3))

    def testPerpendicularTo(self):
        vector = Vector(1, 1, 3)
        vector_z = Vector(0, 0, 1)

        result = vector.perpendicularTo(vector_z)

        self.assertTrue(result == Vector(1, 1, 0))

    def testGetOnePerpendicularVector(self):
        for vector in [Vector(1, 1, 3),
                       Vector(10, 1222, 23),
                       Vector(0.1, 12, -3),
                       Vector(0, 0, 1)]:
            result = vector.getOnePerpendicularVector()
            self.assertAlmostEqual(result.scalarProduct(vector),
                                   0.0)

    def testAngle(self):
        # normalized scalarproduct !!
        vector_x = Vector(1, 0, 0)
        vector_y = Vector(0, 1, 0)
        vector_z = Vector(0, 0, 1)

        vector_xy = Vector(1, 1, 0)

        self.assertAlmostEqual(vector_x.angle(vector_y), pi / 2.0)
        self.assertAlmostEqual(vector_y.angle(vector_z), pi / 2.0)

        self.assertAlmostEqual(vector_x.angle(vector_xy), pi / 4.0)

    def testGetVectorWithAngle(self):

        for vector in [Vector(0, 1, 0),
                       Vector(1, 2, 3),
                       Vector(3, 2, -1)]:
            for angle in arange(0, pi, 0.1):
                vector_with_angle = vector.getVectorWithAngle(angle)
                self.assertAlmostEqual(vector.angle(vector_with_angle),
                                       angle)
