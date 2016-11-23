"""
Unittest for Diffraction class.
"""

import unittest

import numpy as np
import xraylib

from orangecontrib.crystal.diffraction.Diffraction import Diffraction
from orangecontrib.crystal.diffraction.DiffractionSetupSweeps import DiffractionSetupSweeps
from orangecontrib.crystal.diffraction.GeometryType import GeometryType
from orangecontrib.crystal.util.Vector import Vector
from orangecontrib.crystal.util.Photon import Photon
from orangecontrib.crystal.diffraction.GeometryType import BraggDiffraction, LaueDiffraction, BraggTransmission, LaueTransmission
from orangecontrib.crystal.diffraction.DiffractionExceptions import ReflectionImpossibleException, TransmissionImpossibleException, \
                                                                    StructureFactorF0isZeroException, StructureFactorFHisZeroException, \
                                                                    StructureFactorFHbarIsZeroException


class DiffractionTest(unittest.TestCase):

    def assertAlmostEqualLists(self, list1, list2):
        self.assertAlmostEqual(np.linalg.norm(np.array(list1)-np.array(list2)),0,1)

    def assertDiffractionResult(self,energy, s_intensity_fraction, s_phase,p_intensity_fraction, p_phase, diffraction_results):
        self.assertAlmostEqualLists(diffraction_results.sIntensityByEnergy(energy),
                                    s_intensity_fraction)

        self.assertAlmostEqualLists(diffraction_results.sPhaseByEnergy(energy),
                                    s_phase)

        self.assertAlmostEqualLists(diffraction_results.pIntensityByEnergy(energy),
                                    p_intensity_fraction)

        self.assertAlmostEqualLists(diffraction_results.pPhaseByEnergy(energy),
                                    p_phase)

    def testConstructor(self):
        diffraction = Diffraction()
        self.assertIsInstance(diffraction, Diffraction)

    def testCalculateDiffraction(self):
        res = {}
        for geometry_type in GeometryType.allGeometryTypes():
            effective_asymmetry = 0.0
            if geometry_type == LaueDiffraction() or geometry_type == LaueTransmission():
                effective_asymmetry = 90.0

            diffraction_setup = DiffractionSetupSweeps(geometry_type,
                                                 "Si",
                                                 thickness=128 * 1e-6,
                                                 miller_h=1,
                                                 miller_k=1,
                                                 miller_l=1,
                                                 asymmetry_angle=effective_asymmetry,
                                                 energy_min=8174,
                                                 energy_max=8174,
                                                 energy_points=1,
                                                 angle_deviation_min= -20.0e-6,
                                                 angle_deviation_max=20e-6,
                                                 angle_deviation_points=5)
            diffraction = Diffraction()
            res[geometry_type] = diffraction.calculateDiffraction(diffraction_setup)

    def testCalculateBraggDiffraction(self):
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
                                             angle_deviation_min= -20.0e-6,
                                             angle_deviation_max=20e-6,
                                             angle_deviation_points=5)

        diffraction = Diffraction()
        res = diffraction.calculateDiffraction(diffraction_setup)

        s_intensity_fraction=[0.017519141613069177, 0.0321954521714361, 0.07981125895068454, 0.920965084591721, 0.9417181994525138]
        s_phase=[-0.745427562155594, -0.8048350757616735, -0.7441070552657782, -1.0347178161614214, -2.353510138419943]
        p_intensity_fraction=[0.014173087736472335, 0.025303154305706777, 0.06615101317795873, 0.5244213525516417, 0.9369357917670563]
        p_phase=[-0.793312359389805, -0.7582549664194022, -0.750381901971316, -0.8168058020223106, -2.353282699138147]

        self.assertDiffractionResult(10000,
                                     s_intensity_fraction,
                                     s_phase,
                                     p_intensity_fraction,
                                     p_phase,
                                     res)

    def testCalculateBraggTransmission(self):
        diffraction_setup = DiffractionSetupSweeps(BraggTransmission(),
                                             "Si",
                                             thickness=7 * 1e-6,
                                             miller_h=1,
                                             miller_k=1,
                                             miller_l=1,
                                             asymmetry_angle= -5,
                                             energy_min=10174,
                                             energy_max=10174,
                                             energy_points=1,
                                             angle_deviation_min= -20.0e-6,
                                             angle_deviation_max=20e-6,
                                             angle_deviation_points=5)

        diffraction = Diffraction()
        res = diffraction.calculateDiffraction(diffraction_setup)

        s_intensity_fraction=[0.6226567465900791, 0.6438109466925752, 0.6414813069615722, 0.5966674813771604, 0.45178497063185913]
        s_phase=[2.286827125757465, 2.11586718740292, 1.8761281776985377, 1.444935411854202, -0.015769881275207204]
        p_intensity_fraction=[0.6287809489878944, 0.6436830110383608, 0.6260332041734042, 0.5556946212761588, 0.4666570232587092]
        p_phase=[2.4244705128134725, 2.2877506323333496, 2.093850209325308, 1.7465537434885796, 0.8969740263938913]

        self.assertDiffractionResult(10174,
                                     s_intensity_fraction,
                                     s_phase,
                                     p_intensity_fraction,
                                     p_phase,
                                     res)

    def testCalculateLaueDiffraction(self):
        diffraction_setup = DiffractionSetupSweeps(LaueDiffraction(),
                                             "Si",
                                             thickness=100 * 1e-6,
                                             miller_h=1,
                                             miller_k=1,
                                             miller_l=1,
                                             asymmetry_angle=90,
                                             energy_min=8000,
                                             energy_max=8000,
                                             energy_points=1,
                                             angle_deviation_min= -20.0e-6,
                                             angle_deviation_max=20.0e-6,
                                             angle_deviation_points=5)
        diffraction = Diffraction()
        res = diffraction.calculateDiffraction(diffraction_setup)

        s_intensity_fraction=[0.0953161518048925, 0.158471134649239, 0.2844237578381098, 0.158487539849245, 0.09531815291902448]
        s_phase=[2.7878364694515985, -0.816280378494231, -1.6227539168093197, -2.0061870787600458, 0.4081575143878531]
        p_intensity_fraction=[0.0067872399580799405, 0.09329690887082268, 0.12605693490089803, 0.09327296207883676, 0.006786852383095909]
        p_phase=[-1.843856553406182, 1.687240781547736, 0.9198814442403762, 0.49730800506928513, 2.059512850321714]

        self.assertDiffractionResult(8000,
                                     s_intensity_fraction,
                                     s_phase,
                                     p_intensity_fraction,
                                     p_phase,
                                     res)

    def testCalculateLaueTransmission(self):
        diffraction_setup = DiffractionSetupSweeps(LaueTransmission(),
                                             "Si",
                                             thickness=100 * 1e-6,
                                             miller_h=1,
                                             miller_k=1,
                                             miller_l=1,
                                             asymmetry_angle=90,
                                             energy_min=10000,
                                             energy_max=10000,
                                             energy_points=1,
                                             angle_deviation_min= -20.0e-6,
                                             angle_deviation_max=20.0e-6,
                                             angle_deviation_points=5)
        diffraction = Diffraction()
        res = diffraction.calculateDiffraction(diffraction_setup)

        s_intensity_fraction=[0.500009760116572, 0.3730481560147652, 0.1926195176946302, 0.1283757246156211, 0.25695819698222316]
        s_phase=[2.2281966144788545, -0.23994912028908538, -0.215722969718611, 0.25956794505611297, -1.8920272377134075]
        p_intensity_fraction=[0.44963571348593884, 0.5762774883052565, 0.4809772356165785, 0.345952433909957, 0.23751769111657017]
        p_phase=[2.8624375781774436, 0.5308696618055758, 0.1704734474721342, -0.3129214909448153, -2.5856672658533006]

        self.assertDiffractionResult(10000,
                                     s_intensity_fraction,
                                     s_phase,
                                     p_intensity_fraction,
                                     p_phase,
                                     res)

    def testCalculatePsiFromStructureFactor(self):
        diffraction = Diffraction()
        crystal = xraylib.Crystal_GetCrystal("Si")
        photon_in = Photon(8000, Vector(-1,0,-1))
        structure_factor = 113.581288  + 1.763808j

        unitcell_volume = crystal['volume'] * 10 ** -30
        psi = diffraction._calculatePsiFromStructureFactor(unitcell_volume, photon_in, structure_factor)
        self.assertAlmostEqual(psi.real,-1.527826e-5)
        self.assertAlmostEqual(psi.imag,-2.372566e-7)

    def testCheckSetup(self):
        diffraction = Diffraction()

        diffraction_setup = DiffractionSetupSweeps(BraggDiffraction(),
                                             "Si",
                                             thickness=128 * 1e-6,
                                             miller_h=1,
                                             miller_k=1,
                                             miller_l=1,
                                             asymmetry_angle=0,
                                             energy_min=8174 ,
                                             energy_max=8174 ,
                                             energy_points=1 ,
                                             angle_deviation_min= -20.0e-6,
                                             angle_deviation_max=20e-6,
                                             angle_deviation_points=5)

        angle_bragg = 0.19902705045
        F_0     = 113.581288 +  1.763808j
        F_H     =  43.814631 - 42.050823J
        F_H_bar =  42.050823 + 43.814631j

        # Test possible setup.
        diffraction._checkSetup(diffraction_setup,
                                angle_bragg,
                                F_0,
                                F_H,
                                F_H_bar)

        # Test impossible Bragg reflection.
        diffraction_setup._asymmetry_angle = 45

        self.assertRaises(ReflectionImpossibleException,diffraction._checkSetup,diffraction_setup,
                          angle_bragg,F_0,F_H,F_H_bar)

        diffraction_setup._geometry_type = BraggTransmission()
        self.assertRaises(ReflectionImpossibleException,diffraction._checkSetup,diffraction_setup,
                          angle_bragg,F_0,F_H,F_H_bar)

        # Test impossible Laue reflection.
        diffraction_setup._asymmetry_angle = 10

        diffraction_setup._geometry_type = LaueDiffraction()
        self.assertRaises(TransmissionImpossibleException,diffraction._checkSetup,diffraction_setup,
                          angle_bragg,F_0,F_H,F_H_bar)

        diffraction_setup._geometry_type = LaueTransmission()
        self.assertRaises(TransmissionImpossibleException,diffraction._checkSetup,diffraction_setup,
                          angle_bragg,F_0,F_H,F_H_bar)

        # Test forbidden reflection
        diffraction_setup._geometry_type = BraggDiffraction()
        diffraction_setup._asymmetry_angle = 0

        # ... for F_0.
        self.assertRaises(StructureFactorF0isZeroException,diffraction._checkSetup,diffraction_setup,
                          angle_bragg,0.0,F_H,F_H_bar)

        # ... for F_H.
        self.assertRaises(StructureFactorFHisZeroException,diffraction._checkSetup,diffraction_setup,
                          angle_bragg,F_0,0.0,F_H_bar)

        self.assertRaises(StructureFactorFHisZeroException,diffraction._checkSetup,diffraction_setup,
                          angle_bragg,F_0,float('nan')*1j,F_H_bar)

        # ... for F_H_bar.
        self.assertRaises(StructureFactorFHbarIsZeroException,diffraction._checkSetup,diffraction_setup,
                          angle_bragg,F_0,F_H,0.0)

        self.assertRaises(StructureFactorFHbarIsZeroException,diffraction._checkSetup,diffraction_setup,
                          angle_bragg,F_0,F_H,float('nan')*1j)

    @unittest.skip("Do not test against XRT")
    def testXRTDriver(self):
        import orangecontrib.crystal.util.XRTDriver as XRTDriver
        from pylab import plot, legend, ylabel, xlabel, title, savefig, figure

        energy = 8100.0
        for geo in allGeometryTypes():
            for asymmetry in [10.0, 0.0, 5.0]:

                res = XRTDriver.calculateDiffraction(E0=energy,
                                                     alpha=asymmetry)

                xrt_res = res[geo]

                effective_asymmetry = asymmetry
                if geo == LaueDiffraction() or geo == LaueTransmission():
                    effective_asymmetry = 90.0-asymmetry

                diffraction_setup = DiffractionSetupSweeps(geo,
                                                     "Si",
                                                     thickness=100 * 1e-6,
                                                     miller_h=1,
                                                     miller_k=1,
                                                     miller_l=1,
                                                     asymmetry_angle=effective_asymmetry,
                                                     energy=energy,
                                                     angle_deviation_min= -100e-6,
                                                     angle_deviation_max=100e-6,
                                                     angle_deviation_points=300)
                diffraction = Diffraction()
                res = diffraction.calculateDiffraction(diffraction_setup)

                x = [i * 1e+6 for i in res.angleDeviations()]
                plot(x, res.sIntensityByEnergy(), label="S polarization")
                x = [i * 1e+6 for i in xrt_res.angleDeviations()]
                plot(x, xrt_res.sIntensityByEnergy(), label="XRT S polarization")
                legend()
                title(geo.description())
                ylabel('Reflectivity')
                xlabel("Angle deviation in urad")
                filename = "%s_Asym%i_Reflectivity_S.png" % (geo.description().replace(" ", "_"),
                                                             asymmetry)
                savefig(filename)
                figure()

                x = [i * 1e+6 for i in res.angleDeviations()]
                plot(x, res.pIntensityByEnergy(), label="P polarization")
                x = [i * 1e+6 for i in xrt_res.angleDeviations()]
                plot(x, xrt_res.pIntensityByEnergy(), label="XRT P polarization")
                legend()
                title(geo.description())
                ylabel('Reflectivity')
                xlabel("Angle deviation in urad")
                filename = "%s_Asym%i_Reflectivity_P.png" % (geo.description().replace(" ", "_"),
                                                             asymmetry)
                savefig(filename)
                figure()

                x = [i * 1e+6 for i in res.angleDeviations()]
                plot(x, res.sPhaseByEnergy(), label="S polarization")
                x = [i * 1e+6 for i in xrt_res.angleDeviations()]
                plot(x, xrt_res.sPhaseByEnergy(), label="XRT S polarization")
                legend()
                title(geo.description())
                ylabel('Phase shift')
                xlabel("Angle deviation in urad")
                filename = "%s_Asym%i_Phase_S.png" % (geo.description().replace(" ", "_"),
                                                             asymmetry)
                savefig(filename)
                figure()

                x = [i * 1e+6 for i in res.angleDeviations()]
                plot(x, res.pPhaseByEnergy(), label="P polarization")
                x = [i * 1e+6 for i in xrt_res.angleDeviations()]
                plot(x, xrt_res.pPhaseByEnergy(), label="XRT P polarization")
                legend()
                title(geo.description())
                ylabel('Phase shift')
                xlabel("Angle deviation in urad")
                filename = "%s_Asym%i_Phase_P.png" % (geo.description().replace(" ", "_"),
                                                             asymmetry)
                savefig(filename)
                figure()

    @unittest.skip("Do not produce former bug output.")
    def testBugsByLaurence(self):
        geometries = [ BraggTransmission(), LaueTransmission()]
        thicknessses = [128 * 1e-6, 5*1e-6]
        crystal_names = ["Diamond","Si"]
        asymmetries = [0,10,30,50]
        
        plots = []
        for thickness in thicknessses:
            for crystal_name in crystal_names:
                for asymmetry in asymmetries:
                    effective_asymmetry = asymmetry
                    for geo in geometries:
                        if geo == LaueDiffraction() or geo == LaueTransmission():
                            effective_asymmetry = 90.0-asymmetry
                        
                        diffraction_setup = DiffractionSetupSweeps(geo,
                                                     crystal_name,
                                                     thickness=thickness,
                                                     miller_h=1,
                                                     miller_k=1,
                                                     miller_l=1,
                                                     asymmetry_angle=effective_asymmetry,
                                                     energy=3124,
                                                     angle_deviation_min= -120e-6,
                                                     angle_deviation_max=120e-6,
                                                     angle_deviation_points=300)
                        
                        diffraction = Diffraction()
                        try:
                            res = diffraction.calculateDiffraction(diffraction_setup)
                            for p in res.asPlotData1D():
                                plots.append(p)
                        except Exception as ex:
                            print(ex)

        import sys
        from orangecontrib.crystal.widgets.diffraction.PlotViewer1D import PlotViewer1D
        from PyQt4.Qt import QApplication

        application = QApplication(sys.argv)
        ow = PlotViewer1D()
        ow.show()
        ow.setPlots(plots)  
        application.exec_()