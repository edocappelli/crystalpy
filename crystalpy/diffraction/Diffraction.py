"""
Calculates a crystal diffraction.
Except for energy all units are in SI. Energy is in eV.
"""

from math import isnan

from numpy import pi
import scipy.constants.codata

from orangecontrib.crystal.diffraction.GeometryType import BraggDiffraction, BraggTransmission, LaueDiffraction, LaueTransmission
from orangecontrib.crystal.diffraction.DiffractionExceptions import ReflectionImpossibleException, TransmissionImpossibleException, \
                                                                    StructureFactorF0isZeroException, StructureFactorFHisZeroException, \
                                                                    StructureFactorFHbarIsZeroException
from orangecontrib.crystal.util.Photon import Photon
from orangecontrib.crystal.diffraction.DiffractionResult import DiffractionResult
from orangecontrib.crystal.diffraction.PerfectCrystalDiffraction import PerfectCrystalDiffraction


class Diffraction(object):
    
    def __init__(self):
        """
        Constructor.
        """
        # Initialize events to being unhandled.
        self.setOnCalculationStart(None)
        self.setOnProgress(None)
        self.setOnCalculationEnd(None)

    def _calculatePsiFromStructureFactor(self, unit_cell_volume, photon_in, structure_factor):
        """
        Calculates the Psi as defined in Zachariasen [3-95].
        :param unit_cell_volume: Volume of the unit cell.
        :param photon_in: Incoming photon.
        :param structure_factor: Structure factor.
        :return: Psi as defined in Zachariasen [3-95].
        """
        codata = scipy.constants.codata.physical_constants
        classical_electron_radius = codata["classical electron radius"][0]

        psi = (-classical_electron_radius * photon_in.wavelength() ** 2 / (pi * unit_cell_volume)) * structure_factor

        return psi

    @staticmethod
    def log(string):
        """
        Logging function.
        :param string: Message to log.
        """
        print(string)

    def logStructureFactors(self, F_0, F_H, F_H_bar):
        """
        Logs the structure factors.
        :param F_0: Structure factor F_0.
        :param F_H: Structure factor F_H.
        :param F_H_bar: Structure factor F_H_bar.
        """
        self.log("f0: (%f , %f)" % (F_0.real, F_0.imag))
        self.log("fH: (%f , %f)" % (F_H.real, F_H.imag))
        self.log("fHbar: (%f , %f)" % (F_H_bar.real, F_H_bar.imag))

    def setOnCalculationStart(self, on_calculation_start):
        """
        Sets handler for calculation start. The handler is called when the calculation starts.
        :param on_calculation_start: Delegate of calculation start event handler.
        """
        self._on_calculation_start = on_calculation_start

    def _onCalculationStart(self):
        """
        Invokes the calculation start event handler if any is registered.
        """
        if self._on_calculation_start is not None:
            self._on_calculation_start()

    def setOnProgress(self, on_progress):
        """
        Sets handler for calculation progression. The handler is called when the calculation progresses.
        :param on_progress: Delegate of calculation progress event handler.
        """
        self._on_progress = on_progress
        
    def _onProgress(self, index, angle_deviation_points):
        """
        Invokes the calculation progress event handler if any is registered.
        """
        if self._on_progress is not None:
            self._on_progress(index, angle_deviation_points)

    def _onProgressEveryTenPercent(self, index, angle_deviation_points):
        """
        Raises on progress event every ten percent of progression.
        :param index: Index of current point.
        :param angle_deviation_points: Number of total points to calculate.
        """
        ten_percent = angle_deviation_points / 10

        if((index + 1) % ten_percent == 0 or
                index == angle_deviation_points):

            self._onProgress(index + 1, angle_deviation_points)

    def setOnCalculationEnd(self, on_calculation_end):
        """
        Sets handler for calculation end. The handler is called when the calculation ends.
        :param on_calculation_end: Delegate of calculation end event handler.
        """
        self._on_calculation_end = on_calculation_end

    def _onCalculationEnd(self):
        """
        Invokes the calculation end event handler if any is registered.
        """
        if self._on_calculation_end is not None:
            self._on_calculation_end()

    def _checkSetup(self, diffraction_setup, bragg_angle, F_0, F_H, F_H_bar):
        """
        Checks if a given diffraction setup is possible, i.e. if a given Diffraction/Transmission for the given asymmetry
        and Miller indices is possible. Raises an exception if impossible.
        :param diffraction_setup: Diffraction setup.
        :param bragg_angle: Bragg angle.
        :param F_0: Structure factor F_0.
        :param F_H: Structure factor F_H.
        :param F_H_bar: Structure factor F_H_bar.
        """
        # Check if the given geometry is a valid Bragg/Laue geometry.
        if diffraction_setup.geometryType() == BraggDiffraction() or diffraction_setup.geometryType() == BraggTransmission():
            if diffraction_setup.asymmetryAngle() >= bragg_angle:
                raise ReflectionImpossibleException()
        elif diffraction_setup.geometryType() == LaueDiffraction() or diffraction_setup.geometryType() == LaueTransmission():
            if diffraction_setup.asymmetryAngle() <= bragg_angle:
                raise TransmissionImpossibleException()

        # Check structure factor F_0.
        if abs(F_0.real) < 1e-7 or isnan(F_0.real):
            raise StructureFactorF0isZeroException()

        # Check structure factor F_H.
        if abs(F_H.real) < 1e-7 or isnan(F_H.real) or abs(F_H.imag) < 1e-7 or isnan(F_H.imag):
            raise StructureFactorFHisZeroException()

        # Check structure factor F_H_bar.
        if abs(F_H_bar.real) < 1e-7 or isnan(F_H_bar.real) or abs(F_H_bar.imag) < 1e-7 or isnan(F_H_bar.imag):
            raise StructureFactorFHbarIsZeroException()

    def _perfectCrystalForEnergy(self, diffraction_setup, energy):

        # Retrieve bragg angle.
        angle_bragg = diffraction_setup.angleBragg(energy)

        # Get structure factors for all relevant lattice vectors 0,H,H_bar.
        F_0 = diffraction_setup.F0(energy)
        F_H = diffraction_setup.FH(energy)
        F_H_bar = diffraction_setup.FH_bar(energy)

        # Check if given Bragg/Laue geometry and given miller indices are possible.
        self._checkSetup(diffraction_setup, angle_bragg, F_0, F_H, F_H_bar)

        # Log the structure factors.
        self.logStructureFactors(F_0, F_H, F_H_bar)

        # Retrieve lattice spacing d.
        d_spacing = diffraction_setup.dSpacing() * 1e-10

        # Calculate the Bragg normal B_H.
        normal_bragg = diffraction_setup.normalBragg()

        # Calculate the surface normal n.
        normal_surface = diffraction_setup.normalSurface()

        # Calculate the incoming photon direction (parallel to k_0).
        photon_direction = diffraction_setup.incomingPhotonDirection(energy, 0.0)

        # Create photon k_0.
        photon_in = Photon(energy, photon_direction)

        # Retrieve unitcell volume from xraylib.
        unitcell_volume = diffraction_setup.unitcellVolume() * 10 ** -30

        # Calculate psis as defined in Zachariasen [3-95]
        psi_0     = self._calculatePsiFromStructureFactor(unitcell_volume, photon_in, F_0)
        psi_H     = self._calculatePsiFromStructureFactor(unitcell_volume, photon_in, F_H)
        psi_H_bar = self._calculatePsiFromStructureFactor(unitcell_volume, photon_in, F_H_bar)

        # Create PerfectCrystalDiffraction instance.
        perfect_crystal = PerfectCrystalDiffraction(geometry_type=diffraction_setup.geometryType(),
                                                    bragg_normal=normal_bragg,
                                                    surface_normal=normal_surface,
                                                    bragg_angle=angle_bragg,
                                                    psi_0=psi_0,
                                                    psi_H=psi_H,
                                                    psi_H_bar=psi_H_bar,
                                                    thickness=diffraction_setup.thickness(),
                                                    d_spacing=d_spacing)

        return perfect_crystal

    def _calculateDiffractionForEnergy(self, diffraction_setup, energy, result):
        """
        Calculates the diffraction/transmission given by the setup.
        :param diffraction_setup: The diffraction setup.
        :return: DiffractionResult representing this setup.
        """
        # Get PerfectCrystal instance for the current energy.
        perfect_crystal = self._perfectCrystalForEnergy(diffraction_setup, energy)

        # Raise calculation start.
        self._onCalculationStart()

        # For every deviation from Bragg angle ...
        for index, deviation in enumerate(diffraction_setup.angleDeviationGrid()):

            # Raise OnProgress event if progressed by 10 percent.
            self._onProgressEveryTenPercent(index, diffraction_setup.angleDeviationPoints())

            # Calculate deviated incoming photon.
            photon_direction = diffraction_setup.incomingPhotonDirection(energy, deviation)
            photon_in = Photon(energy, photon_direction)

            # Calculate diffraction for current incoming photon.
            result_deviation = perfect_crystal.calculateDiffraction(photon_in)

            # Calculate polarization difference between sigma and pi polarization.
            polarization_difference = result_deviation["S"] / result_deviation["P"]

            # Add result of current deviation.
            result.add(energy,
                       deviation,
                       result_deviation["S"],
                       result_deviation["P"],
                       polarization_difference)

        # Raise calculation end.
        self._onCalculationEnd()

        # Return diffraction results.
        return result

    def calculateDiffraction(self, diffraction_setup):
        """
        Calculates the diffraction/transmission given by the setup.
        :param diffraction_setup: The diffraction setup.
        :return: DiffractionResult representing this setup.
        """
        # Create DiffractionResult instance.
        result = DiffractionResult(diffraction_setup, 0.0)

        for energy in diffraction_setup.energies():
            self._calculateDiffractionForEnergy(diffraction_setup, energy, result)

        # Return diffraction results.
        return result
