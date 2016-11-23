"""
Unittest for StokesVector class.
"""
import numpy as np

import unittest

from orangecontrib.crystal.polarization.StokesVector import StokesVector


def _generate_element_list():
    element_list = [0.78177969457877930,
                    0.22595711869558588,
                    0.28797567756487550,
                    0.58551861060989900]
    return element_list


class StokesVectorTest(unittest.TestCase):
    def setUp(self):
        self.element_list = _generate_element_list()
        self.stokes_vector = StokesVector(self.element_list)

    def test_constructor(self):
        self.assertIsInstance(self.stokes_vector, StokesVector)
        self.assertEqual(self.stokes_vector.s0, 0.78177969457877930)
        self.assertEqual(self.stokes_vector.s1, 0.22595711869558588)
        self.assertEqual(self.stokes_vector.s2, 0.28797567756487550)
        self.assertEqual(self.stokes_vector.s3, 0.58551861060989900)

    def test_get_array(self):
        array1 = self.stokes_vector.get_array(numpy=True)
        array2 = self.stokes_vector.get_array(numpy=False)

        self.assertEqual(type(array1), np.ndarray)
        self.assertEqual(type(array2), list)
        np.testing.assert_array_equal(array1, np.asarray(self.element_list))
        self.assertListEqual(array2, self.element_list)

    def test_polarization_degree(self):
        pol_deg = self.stokes_vector.polarization_degree()
        self.assertEqual(type(pol_deg), float)
        self.assertAlmostEqual(pol_deg, 0.7489560226111716)

    def test_operator_equal(self):
        stokes_vector1 = self.stokes_vector
        stokes_vector2 = StokesVector([0.7817796945787793,
                                       0.22595711869558588,
                                       0.2879756775648755,
                                       0.585518610609899])  # without final zeros
        stokes_vector3 = StokesVector([round(0.78177969457877930, 6),
                                       round(0.22595711869558588, 6),
                                       round(0.28797567756487550, 6),
                                       round(0.58551861060989900, 6)])  # rounded float
        self.assertTrue(stokes_vector1 == stokes_vector1)  # identity
        self.assertTrue(stokes_vector1 == stokes_vector2)
        self.assertFalse(stokes_vector1 == stokes_vector3)

    def test_operator_not_equal(self):
        stokes_vector1 = self.stokes_vector
        stokes_vector2 = StokesVector([0.7817796945787793,
                                       0.22595711869558588,
                                       0.2879756775648755,
                                       0.585518610609899])  # without final zeros
        stokes_vector3 = StokesVector([round(0.78177969457877930, 6),
                                       round(0.22595711869558588, 6),
                                       round(0.28797567756487550, 6),
                                       round(0.58551861060989900, 6)])  # rounded float
        self.assertFalse(stokes_vector1 != stokes_vector1)
        self.assertFalse(stokes_vector1 != stokes_vector2)
        self.assertTrue(stokes_vector1 != stokes_vector3)
