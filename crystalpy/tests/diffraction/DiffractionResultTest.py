"""
Unittest for DiffractionResult class.
"""

import unittest

import numpy

from orangecontrib.crystal.diffraction.ComplexAmplitude import ComplexAmplitude
from orangecontrib.crystal.diffraction.DiffractionSetupSweeps import DiffractionSetupSweeps
from orangecontrib.crystal.diffraction.Diffraction import Diffraction
from orangecontrib.crystal.diffraction.DiffractionResult import DiffractionResult
from orangecontrib.crystal.diffraction.GeometryType import BraggDiffraction


def diffractionSetupSingleEnergy():
    diffraction_setup = DiffractionSetupSweeps(BraggDiffraction(),
                                         "Si",
                                         thickness=0.0100 * 1e-2,
                                         miller_h=1,
                                         miller_k=1,
                                         miller_l=1,
                                         asymmetry_angle=3,
                                         energy_min=10000,
                                         energy_max=10000,
                                         energy_points=1,
                                         angle_deviation_min=-100.0e-6,
                                         angle_deviation_max=100e-6,
                                         angle_deviation_points=50)
    return diffraction_setup


def diffractionSetupMultipleEnergy():
    diffraction_setup = DiffractionSetupSweeps(BraggDiffraction(),
                                         "Si",
                                         thickness=0.0100 * 1e-2,
                                         miller_h=1,
                                         miller_k=1,
                                         miller_l=1,
                                         asymmetry_angle=3,
                                         energy_min=10000,
                                         energy_max=10001,
                                         energy_points=2,
                                         angle_deviation_min=-100.0e-6,
                                         angle_deviation_max=100e-6,
                                         angle_deviation_points=50)
    return diffraction_setup


def diffractionResult():
    diffraction_setup = diffractionSetupSingleEnergy()

    bragg_angle = 0.3
    diffraction_result = DiffractionResult(diffraction_setup,
                                           bragg_angle)
    return diffraction_result


class DiffractionResultTest(unittest.TestCase):
    def testConstructor(self):
        diffraction_setup = diffractionSetupSingleEnergy()

        bragg_angle = 0.3
        diffraction_result = DiffractionResult(diffraction_setup,
                                               bragg_angle)

        self.assertIsInstance(diffraction_result,
                              DiffractionResult)

        self.assertEqual(diffraction_result._diffraction_setup,
                         diffraction_setup)

        self.assertEqual(diffraction_result._bragg_angle,
                         bragg_angle)

        self.assertEqual(diffraction_result.angleDeviations().shape[0],
                         50)

        self.assertEqual(diffraction_result._intensities.shape,
                         (1, 50, 3))
        self.assertEqual(diffraction_result._phases.shape,
                         (1, 50, 3))

    def testDiffractionSetup(self):
        diffraction_setup = diffractionSetupSingleEnergy()

        bragg_angle = 0.3
        diffraction_result = DiffractionResult(diffraction_setup,
                                               bragg_angle)

        self.assertEqual(diffraction_result.diffractionSetup(),
                         diffraction_setup)


    def testBraggAngle(self):
        diffraction_setup = diffractionSetupSingleEnergy()

        bragg_angle = 0.3
        diffraction_result = DiffractionResult(diffraction_setup,
                                               bragg_angle)
        self.assertEqual(diffraction_result.braggAngle(),
                         bragg_angle)

    def testAngleDeviations(self):
        diffraction_result = diffractionResult()
        self.assertAlmostEqual(diffraction_result.angleDeviations()[0],-100e-6)
        self.assertAlmostEqual(diffraction_result.angleDeviations()[24],-2.040816e-06)
        self.assertAlmostEqual(diffraction_result.angleDeviations()[25],2.040816e-06)
        self.assertAlmostEqual(diffraction_result.angleDeviations()[49],100e-6)

    def testAngles(self):
        diffraction_result = diffractionResult()
        self.assertAlmostEqual(diffraction_result.angles()[0], 0.3+-100e-6)
        self.assertAlmostEqual(diffraction_result.angles()[24], 0.3+-2.040816e-06)
        self.assertAlmostEqual(diffraction_result.angles()[25], 0.3+2.040816e-06)
        self.assertAlmostEqual(diffraction_result.angles()[49], 0.3+100e-6)

    def testSIntensityByEnergy(self):
        diffraction_result = diffractionResult()
        self.assertTrue((diffraction_result.sIntensityByEnergy(10000)==numpy.zeros(50)).all())

    def testSPhaseByEnergy(self):
        diffraction_result = diffractionResult()
        self.assertTrue((diffraction_result.sPhaseByEnergy(10000)==numpy.zeros(50)).all())

    def testPIntensityByEnergy(self):
        diffraction_result = diffractionResult()
        self.assertTrue((diffraction_result.pIntensityByEnergy(10000)==numpy.zeros(50)).all())

    def testPPhaseByEnergy(self):
        diffraction_result = diffractionResult()
        self.assertTrue((diffraction_result.pPhaseByEnergy(10000)==numpy.zeros(50)).all())

    def testDifferenceIntensityByEnergy(self):
        diffraction_result = diffractionResult()
        self.assertTrue((diffraction_result.differenceIntensityByEnergy(10000)==numpy.zeros(50)).all())

    def testDifferencePhaseByEnergy(self):
        diffraction_result = diffractionResult()
        self.assertTrue((diffraction_result.differencePhaseByEnergy(10000)==numpy.zeros(50)).all())

    def testSIntensityByDeviation(self):
        diffraction_result = diffractionResult()
        self.assertTrue((diffraction_result.sIntensityByDeviation(100e6)==numpy.zeros(1)).all())

    def testSPhaseByDeviation(self):
        diffraction_result = diffractionResult()
        self.assertTrue((diffraction_result.sPhaseByDeviation(100e6)==numpy.zeros(1)).all())

    def testPIntensityByDeviation(self):
        diffraction_result = diffractionResult()
        self.assertTrue((diffraction_result.pIntensityByDeviation(100e6)==numpy.zeros(1)).all())

    def testPPhaseByDeviation(self):
        diffraction_result = diffractionResult()
        self.assertTrue((diffraction_result.pPhaseByDeviation(100e6)==numpy.zeros(1)).all())

    def testDifferenceIntensityByDeviation(self):
        diffraction_result = diffractionResult()
        self.assertTrue((diffraction_result.differenceIntensityByDeviation(100e6)==numpy.zeros(1)).all())

    def testDifferencePhaseByDeviation(self):
        diffraction_result = diffractionResult()
        self.assertTrue((diffraction_result.differencePhaseByDeviation(100e6)==numpy.zeros(1)).all())

    def testAdd(self):
        diffraction_setup = diffractionSetupMultipleEnergy()

        bragg_angle = 0.3
        diffraction_result = DiffractionResult(diffraction_setup,
                                               bragg_angle)

        for energy in (10000, 10001):
            for deviation in (-100e-6, -2.040816e-06, 2.040816e-06, 100e-6):
                s_complex_amplitude = ComplexAmplitude(energy+deviation * 1e+4)
                p_complex_amplitude = ComplexAmplitude(energy+2.0*deviation * 1e+4)
                difference_complex_amplitude = s_complex_amplitude / p_complex_amplitude

                diffraction_result.add(energy,
                                       deviation,
                                       s_complex_amplitude,
                                       p_complex_amplitude,
                                       difference_complex_amplitude)

        # Test for energy 10000
        # s polarization
        self.assertAlmostEqual(diffraction_result.pIntensityByEnergy(10000)[0],
                               99960004.0)
        self.assertAlmostEqual(diffraction_result.pIntensityByEnergy(10000)[24],
                               99999183.675265953)
        self.assertAlmostEqual(diffraction_result.pIntensityByEnergy(10000)[25],
                               100000816.32806599)
        self.assertAlmostEqual(diffraction_result.pIntensityByEnergy(10000)[49],
                               100040004.0)

        # p polarization
        self.assertAlmostEqual(diffraction_result.sIntensityByEnergy(10000)[0],
                               99980001.0)
        self.assertAlmostEqual(diffraction_result.sIntensityByEnergy(10000)[24],
                               99999591.837216482)
        self.assertAlmostEqual(diffraction_result.sIntensityByEnergy(10000)[25],
                               100000408.16361651)
        self.assertAlmostEqual(diffraction_result.sIntensityByEnergy(10000)[49],
                               100020001.0)


        # Test for energy 10001
        # s polarization
        self.assertAlmostEqual(diffraction_result.sIntensityByEnergy(10001)[0],
                               100000000.0)
        self.assertAlmostEqual(diffraction_result.sIntensityByEnergy(10001)[24],
                               100019592.79640016)
        self.assertAlmostEqual(diffraction_result.sIntensityByEnergy(10001)[25],
                               100020409.20443282)
        self.assertAlmostEqual(diffraction_result.sIntensityByEnergy(10001)[49],
                               100040004.0)

        # p polarization
        self.assertAlmostEqual(diffraction_result.pIntensityByEnergy(10001)[0],
                               99980001.0)
        self.assertAlmostEqual(diffraction_result.pIntensityByEnergy(10001)[24],
                               100019184.59363331)
        self.assertAlmostEqual(diffraction_result.pIntensityByEnergy(10001)[25],
                               100020817.40969864)
        self.assertAlmostEqual(diffraction_result.pIntensityByEnergy(10001)[49],
                               100060009.0)

    def testAsPlotData1D(self):
        diffraction_setup = DiffractionSetupSweeps(BraggDiffraction(),
                                             "Si",
                                             thickness=0.0100 * 1e-2,
                                             miller_h=1,
                                             miller_k=1,
                                             miller_l=1,
                                             asymmetry_angle=3,
                                             energy_min=10000,
                                             energy_max=10000,
                                             energy_points=1,
                                             angle_deviation_min= -100.0e-6,
                                             angle_deviation_max=100e-6,
                                             angle_deviation_points=50)

        diffraction = Diffraction()
        res = diffraction.calculateDiffraction(diffraction_setup)

        plot_data_2d = res.asPlotData1D()
        plot_info =  plot_data_2d[0].plotInfo()
        self.assertIn("Geometry Type", plot_info)
        self.assertIn("Crystal Name", plot_info)
        self.assertIn("Miller indices (h,k,l)", plot_info)

        self.assertEqual(plot_info["Geometry Type"], "Bragg diffraction")
        self.assertEqual(plot_info["Crystal Name"], "Si")
        self.assertEqual(plot_info["Miller indices (h,k,l)"], "(1,1,1)")