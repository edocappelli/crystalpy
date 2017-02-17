"""
Unittest for MuellerMatrix class.
"""
import unittest

import numpy as np

from numpy.testing import assert_array_almost_equal

from crystalpy.polarization.MuellerMatrix import MuellerMatrix


def _generate_matrix():
    matrix = np.asarray([[0.3259811941363535, 0.47517601596586845, 0.5518825026595213, 0.309697589767706],
                         [0.9302512754240655, 0.31523109368760793, 0.14892732688907617, 0.6794292697143903],
                         [0.4445952221208559, 0.3645127814592628, 0.1683036723562683, 0.6726938039088263],
                         [0.5613952867979658, 0.1067151292280123, 0.4887844163526709, 0.27447059819338904]])
    return matrix


class MuellerMatrixTest(unittest.TestCase):
    def setUp(self):
        self.matrix = _generate_matrix()
        self.mueller_matrix = MuellerMatrix(self.matrix)

    def test_constructor(self):
        matrix = self.mueller_matrix.matrix
        self.assertEqual(type(matrix), np.ndarray)
        self.assertListEqual(matrix.flatten().tolist(), self.matrix.flatten().tolist())

    def test_from_matrix_to_elements(self):
        element_list = [0.3259811941363535, 0.47517601596586845, 0.5518825026595213, 0.309697589767706,
                        0.9302512754240655, 0.31523109368760793, 0.14892732688907617, 0.6794292697143903,
                        0.4445952221208559, 0.3645127814592628, 0.1683036723562683, 0.6726938039088263,
                        0.5613952867979658, 0.1067151292280123, 0.4887844163526709, 0.27447059819338904]
        matrix1 = self.mueller_matrix.from_matrix_to_elements(return_numpy=False)
        matrix2 = self.mueller_matrix.from_matrix_to_elements(return_numpy=True)
        self.assertEqual(type(matrix1), list)
        self.assertEqual(type(matrix2), np.ndarray)
        np.testing.assert_array_almost_equal(matrix2, np.asarray(element_list))
        self.assertListEqual(matrix1, element_list)

    def test_matrix_by_scalar(self):
        matrix = self.mueller_matrix.matrix
        scalar1 = 2.2562276967724735
        scalar2 = scalar1 * 1e-8
        scalar3 = -5.010632784090829
        self.assertIsInstance(self.mueller_matrix.matrix_by_scalar(scalar1), MuellerMatrix)
        # once I have checked the products are MuellerMatrix objects, I use the __eq__ method.
        self.assertTrue(self.mueller_matrix.matrix_by_scalar(scalar1) == MuellerMatrix(matrix * scalar1))
        self.assertIsInstance(self.mueller_matrix.matrix_by_scalar(scalar2), MuellerMatrix)
        self.assertTrue(self.mueller_matrix.matrix_by_scalar(scalar2) == MuellerMatrix(matrix * scalar2))
        self.assertIsInstance(self.mueller_matrix.matrix_by_scalar(scalar3), MuellerMatrix)
        self.assertTrue(self.mueller_matrix.matrix_by_scalar(scalar3) == MuellerMatrix(matrix * scalar3))

    def test_matrix_by_vector(self):
        vector = np.array([2.3305949456965642, 8.395751778702131, 1.8896988291928611, -9.515370772747518])
        res1 = self.mueller_matrix.matrix_by_vector(vector, return_numpy=False)
        res2 = self.mueller_matrix.matrix_by_vector(vector, return_numpy=True)
        self.assertEqual(type(res1), list)
        self.assertEqual(type(res2), np.ndarray)
        self.assertListEqual(res1, list(np.dot(self.mueller_matrix.matrix, vector)))
        np.testing.assert_array_almost_equal(res2, np.dot(self.mueller_matrix.matrix, vector))

    def test_vector_by_matrix(self):
        vector = np.array([2.3305949456965642, 8.395751778702131, 1.8896988291928611, -9.515370772747518])
        res1 = self.mueller_matrix.vector_by_matrix(vector, return_numpy=False)
        res2 = self.mueller_matrix.vector_by_matrix(vector, return_numpy=True)
        self.assertEqual(type(res1), list)
        self.assertEqual(type(res2), np.ndarray)
        self.assertListEqual(res1, list(np.dot(vector, self.mueller_matrix.matrix)))
        np.testing.assert_array_almost_equal(res2, np.dot(vector, self.mueller_matrix.matrix))

    def test_mueller_times_mueller(self):
        mueller_matrix1 = self.mueller_matrix
        matrix2 = np.asarray([[-6.17688112,  5.02032824,  6.23225893, -0.32032499],
                              [7.27366317,  1.43011512,  1.81816476,  9.62888406],
                              [5.27430283, -1.37894077, -0.99747857,  8.80405558],
                              [1.89547814, -1.17966381, -7.85411672, -7.73368115]])
        product1 = mueller_matrix1.mueller_times_mueller(matrix2, mod=False)
        product2 = mueller_matrix1.mueller_times_mueller(matrix2, mod=True)
        self.assertIsInstance(product1, MuellerMatrix)
        self.assertIsInstance(product2, MuellerMatrix)
        # once I have checked the products are MuellerMatrix objects, I use the __eq__ method.
        self.assertTrue(product1 == MuellerMatrix(np.dot(matrix2, mueller_matrix1.matrix)))
        self.assertTrue(product2 == MuellerMatrix(np.dot(mueller_matrix1.matrix, matrix2)))

    def test_operator_equal(self):
        candidate1 = np.asarray([[0.3259811941363535, 0.47517601596586845, 0.55188250265952130, 0.30969758976770600],
                                [0.9302512754240655, 0.31523109368760793, 0.14892732688907617, 0.67942926971439030],
                                [0.4445952221208559, 0.36451278145926280, 0.16830367235626830, 0.67269380390882630],
                                [0.5613952867979658, 0.10671512922801230, 0.48878441635267090, 0.27447059819338904]])
        candidate1 = MuellerMatrix(candidate1)
        candidate2 = np.asarray([[-6.17688112,  5.02032824,  6.23225893, -0.32032499],
                                [7.27366317,  1.43011512,  1.81816476,  9.62888406],
                                [5.27430283, -1.37894077, -0.99747857,  8.80405558],
                                [1.89547814, -1.17966381, -7.85411672, -7.73368115]])
        candidate2 = MuellerMatrix(candidate2)
        self.assertTrue(self.mueller_matrix == self.mueller_matrix)  # identity
        self.assertTrue(candidate1 == self.mueller_matrix)
        self.assertFalse(candidate2 == self.mueller_matrix)

    def test_operator_not_equal(self):
        candidate1 = np.asarray([[0.3259811941363535, 0.47517601596586845, 0.55188250265952130, 0.30969758976770600],
                                 [0.9302512754240655, 0.31523109368760793, 0.14892732688907617, 0.67942926971439030],
                                 [0.4445952221208559, 0.36451278145926280, 0.16830367235626830, 0.67269380390882630],
                                 [0.5613952867979658, 0.10671512922801230, 0.48878441635267090, 0.27447059819338904]])
        candidate1 = MuellerMatrix(candidate1)
        candidate2 = np.asarray([[-6.17688112, 5.02032824, 6.23225893, -0.32032499],
                                 [7.27366317, 1.43011512, 1.81816476, 9.62888406],
                                 [5.27430283, -1.37894077, -0.99747857, 8.80405558],
                                 [1.89547814, -1.17966381, -7.85411672, -7.73368115]])
        candidate2 = MuellerMatrix(candidate2)
        self.assertFalse(self.mueller_matrix != self.mueller_matrix)  # identity
        self.assertFalse(candidate1 != self.mueller_matrix)
        self.assertTrue(candidate2 != self.mueller_matrix)

    def test_linear_polarizers(self):

        #
        assert_array_almost_equal( MuellerMatrix.initialize_as_linear_polarizer_horizontal().get_matrix(),
                                   0.5 * np.array([
                                                [1.0,1.0,0.0,0.0],
                                                [1.0,1.0,0.0,0.0],
                                                [0.0,0.0,0.0,0.0],
                                                [0.0,0.0,0.0,0.0] ]))

        #
        assert_array_almost_equal( MuellerMatrix.initialize_as_linear_polarizer_horizontal().get_matrix(),
                                   0.5 * np.array([
                                                [1.0,1.0,0.0,0.0],
                                                [1.0,1.0,0.0,0.0],
                                                [0.0,0.0,0.0,0.0],
                                                [0.0,0.0,0.0,0.0] ]))

        #
        assert_array_almost_equal( MuellerMatrix.initialize_as_linear_polarizer_vertical().get_matrix(),
                                   0.5 * np.array([
                                                [1.0,-1.0,0.0,0.0],
                                                [-1.0,1.0,0.0,0.0],
                                                [0.0,0.0,0.0,0.0],
                                                [0.0,0.0,0.0,0.0] ]))


        #
        assert_array_almost_equal( MuellerMatrix.initialize_as_linear_polarizer_plus45().get_matrix(),
                                   0.5 * np.array([
                                                [1.0,0.0,1.0,0.0],
                                                [0.0,0.0,0.0,0.0],
                                                [1.0,0.0,1.0,0.0],
                                                [0.0,0.0,0.0,0.0] ]))


        #
        assert_array_almost_equal( MuellerMatrix.initialize_as_linear_polarizer_minus45().get_matrix(),
                                   0.5 * np.array([
                                                [1.0,0.0,-1.0,0.0],
                                                [0.0,0.0,0.0,0.0],
                                                [-1.0,0.0,1.0,0.0],
                                                [0.0,0.0,0.0,0.0] ]))

    def test_linear_retarders(self):
        #
        assert_array_almost_equal( MuellerMatrix.initialize_as_quarter_wave_plate_fast_vertical().get_matrix(),
                                   np.array([
                                                [1.0,0.0,0.0,0.0],
                                                [0.0,1.0,0.0,0.0],
                                                [0.0,0.0,0.0,-1.0],
                                                [0.0,0.0,1.0,0.0] ]))
        #
        assert_array_almost_equal( MuellerMatrix.initialize_as_quarter_wave_plate_fast_horizontal().get_matrix(),
                                   np.array([
                                                [1.0,0.0,0.0,0.0],
                                                [0.0,1.0,0.0,0.0],
                                                [0.0,0.0,0.0,1.0],
                                                [0.0,0.0,-1.0,0.0] ]))
        #
        assert_array_almost_equal( MuellerMatrix.initialize_as_half_wave_plate().get_matrix(),
                                   np.array([
                                                [1.0,0.0,0.0,0.0],
                                                [0.0,1.0,0.0,0.0],
                                                [0.0,0.0,-1.0,0.0],
                                                [0.0,0.0,0.0,-1.0] ]))
    def test_filter(self):
        #
        assert_array_almost_equal( MuellerMatrix.initialize_as_filter(55.5).get_matrix(),
                                   55.5*np.array([
                                                [1.0,0.0,0.0,0.0],
                                                [0.0,1.0,0.0,0.0],
                                                [0.0,0.0,1.0,0.0],
                                                [0.0,0.0,0.0,1.0] ]))

