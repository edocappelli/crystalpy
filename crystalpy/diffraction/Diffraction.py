"""
Calculates a crystal diffraction.
Except for energy all units are in SI. Energy is in eV.
"""

from math import isnan

from numpy import pi
import scipy.constants.codata

from crystalpy.diffraction.GeometryType import BraggDiffraction, BraggTransmission, LaueDiffraction, LaueTransmission
from crystalpy.diffraction.DiffractionExceptions import ReflectionImpossibleException, TransmissionImpossibleException, \
                                                                    StructureFactorF0isZeroException, StructureFactorFHisZeroException, \
                                                                    StructureFactorFHbarIsZeroException
from crystalpy.util.Photon import Photon
from crystalpy.util.PolarizedPhoton import PolarizedPhoton
from crystalpy.util.ComplexAmplitudePhotonBunch import ComplexAmplitudePhotonBunch

from crystalpy.util.PhotonBunch import PhotonBunch
from crystalpy.util.PolarizedPhotonBunch import PolarizedPhotonBunch


from crystalpy.diffraction.DiffractionResult import DiffractionResult
from crystalpy.diffraction.PerfectCrystalDiffraction import PerfectCrystalDiffraction
from crystalpy.polarization.CrystalPhasePlate import CrystalPhasePlate

from crystalpy.diffraction.DiffractionSetupSweeps import DiffractionSetupSweeps

class Diffraction(object):
    isDebug = False
    
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
        if self.isDebug:
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
        if not isinstance(diffraction_setup,DiffractionSetupSweeps):
            raise Exception("Inmut must be of type: DiffractionSetupSweeps")
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

            # Calculate polarization difference between pi and sigma polarization.
            polarization_difference = result_deviation["P"] / result_deviation["S"]

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

        if not isinstance(diffraction_setup,DiffractionSetupSweeps):
            raise Exception("Input object must be of type DiffractionSetupSweeps")

        # Create DiffractionResult instance.
        result = DiffractionResult(diffraction_setup, 0.0)

        for energy in diffraction_setup.energies():
            self._calculateDiffractionForEnergy(diffraction_setup, energy, result)

        # Return diffraction results.
        return result

    #TODO remove? where is used??

    # def calculateDiffractionForPhotonBunch(self, diffraction_setup, photon_bunch):
    #     """
    #     Calculates the diffraction/transmission given by the setup.
    #     :param diffraction_setup: The diffraction setup.
    #     :return: DiffractionResult representing this setup.
    #     """
    #     if not isinstance(photon_bunch,PhotonBunch):
    #         raise Exception("Input object must be of type PhotonBunch")
    #
    #     # Create DiffractionResult instance.
    #     result = DiffractionResult(diffraction_setup, 0.0)
    #
    #     for energy in diffraction_setup.energies():
    #         self._calculateDiffractionForEnergy(diffraction_setup, energy, result)
    #
    #     # Return diffraction results.
    #     return result

##################################################################################################
# FUNCTIONS ADAPTED TO WORK WITH GENERAL BUNCHES OF PHOTONS AND NOT WITH DIRECTION/ENERGY SWEEPS #
##################################################################################################

    def _perfectCrystalForPhoton(self, diffraction_setup, polarized_photon):

        energy = polarized_photon.energy()

        # Retrieve bragg angle.
        angle_bragg = diffraction_setup.angleBragg(energy)

        # Get structure factors for all relevant lattice vectors 0,H,H_bar.
        F_0 = diffraction_setup.F0(energy)
        F_H = diffraction_setup.FH(energy)
        F_H_bar = diffraction_setup.FH_bar(energy)

        # Check if given Bragg/Laue geometry and given miller indices are possible.
        self._checkSetup(diffraction_setup, angle_bragg, F_0, F_H, F_H_bar)

        # Log the structure factors.
        if self.isDebug:
            self.logStructureFactors(F_0, F_H, F_H_bar)

        # Retrieve lattice spacing d.
        d_spacing = diffraction_setup.dSpacing() * 1e-10

        # Calculate the Bragg normal B_H.
        normal_bragg = diffraction_setup.normalBragg()

        # Calculate the surface normal n.
        normal_surface = diffraction_setup.normalSurface()

        # Retrieve unitcell volume from xraylib.
        unitcell_volume = diffraction_setup.unitcellVolume() * 10 ** -30

        # Calculate psis as defined in Zachariasen [3-95]
        psi_0 = self._calculatePsiFromStructureFactor(unitcell_volume, polarized_photon, F_0)
        psi_H = self._calculatePsiFromStructureFactor(unitcell_volume, polarized_photon, F_H)
        psi_H_bar = self._calculatePsiFromStructureFactor(unitcell_volume, polarized_photon, F_H_bar)

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

    def calculateDiffractedPolarizedPhoton(self, diffraction_setup, incoming_polarized_photon, inclination_angle):
        """
        Calculates the diffraction/transmission given by the setup.
        :param diffraction_setup: The diffraction setup.
        :return: PhotonBunch object made up of diffracted/transmitted photons.
        """
        # Retrieve the incoming Stokes vector.
        incoming_stokes_vector = incoming_polarized_photon.stokesVector()

        # Get PerfectCrystal instance for the current photon.
        perfect_crystal = self._perfectCrystalForPhoton(diffraction_setup, incoming_polarized_photon)

        # Calculate diffraction for current incoming photon.
        complex_amplitudes = perfect_crystal.calculateDiffraction(incoming_polarized_photon)

        # Calculate outgoing Photon.
        outgoing_photon = perfect_crystal._calculatePhotonOut(incoming_polarized_photon)

        # Calculate intensities and phases of the crystal  reflectivities or transmitivities
        intensity_pi = complex_amplitudes["P"].intensity()
        intensity_sigma = complex_amplitudes["S"].intensity()
        phase_pi = complex_amplitudes["P"].phase()
        phase_sigma = complex_amplitudes["S"].phase()

        # Get a CrystalPhasePlate instance which contains the Mueller matrix
        phase_plate = CrystalPhasePlate( #incoming_stokes_vector=incoming_stokes_vector,
                                        intensity_sigma=intensity_sigma,
                                        phase_sigma=phase_sigma,
                                        intensity_pi=intensity_pi,
                                        phase_pi=phase_pi,
                                        inclination_angle=inclination_angle)

        # Use intensities and phases to calculate the Stokes vector for the outgoing photon.
        outgoing_stokes_vector = phase_plate.calculate_stokes_vector(incoming_stokes_vector)

        # Piece together the PolarizedPhoton object.
        outgoing_polarized_photon = PolarizedPhoton(energy_in_ev=outgoing_photon.energy(),
                                                    direction_vector=outgoing_photon.unitDirectionVector(),
                                                    stokes_vector=outgoing_stokes_vector)

        return outgoing_polarized_photon


    def calculateDiffractedPolarizedPhotonBunch(self, diffraction_setup, incoming_bunch, inclination_angle):
        """
        Calculates the diffraction/transmission given by the setup.
        :param diffraction_setup: The diffraction setup.
        :return: PhotonBunch object made up of diffracted/transmitted photons.
        """
        # Create PhotonBunch instance.
        outgoing_bunch = PolarizedPhotonBunch([])

        # Retrieve the photon bunch from the diffraction setup.
        # incoming_bunch = diffraction_setup.incomingPhotons()

        # Check that photon_bunch is indeed a PhotonBunch object.
        if not isinstance(incoming_bunch, PolarizedPhotonBunch):
            raise Exception("The incoming photon bunch must be a PolarizedPhotonBunch object!")

        # Raise calculation start.
        self._onCalculationStart()

        for index, polarized_photon in enumerate(incoming_bunch):

            # Raise OnProgress event if progressed by 10 percent.
            self._onProgressEveryTenPercent(index, len(incoming_bunch))

            outgoing_polarized_photon = self.calculateDiffractedPolarizedPhoton(diffraction_setup, polarized_photon, inclination_angle)
            # Add result of current deviation.
            outgoing_bunch.addPhoton(outgoing_polarized_photon)

        # Raise calculation end.
        self._onCalculationEnd()

        # Return diffraction results.
        return outgoing_bunch

    # calculate complex reflectivity and transmitivity
    def calculateDiffractedComplexAmplitudes(self, diffraction_setup, incoming_photon):

        # Get PerfectCrystal instance for the current photon.
        perfect_crystal = self._perfectCrystalForPhoton(diffraction_setup, incoming_photon)

        # Calculate diffraction for current incoming photon.
        complex_amplitudes = perfect_crystal.calculateDiffraction(incoming_photon)

        return complex_amplitudes

    # using ComplexAmplitudePhoton
    def calculateDiffractedComplexAmplitudePhoton(self, diffraction_setup,photon):

        # Get PerfectCrystal instance for the current photon.
        perfect_crystal = self._perfectCrystalForPhoton(diffraction_setup, photon)

        coeffs = self.calculateDiffractedComplexAmplitudes(diffraction_setup,photon)

        # Calculate outgoing Photon.
        outgoing_photon = perfect_crystal._calculatePhotonOut(photon)
        # apply reflectivities
        outgoing_photon.rescaleEsigma(coeffs["S"])
        outgoing_photon.rescaleEpi(coeffs["P"])

        return outgoing_photon

    def calculateDiffractedComplexAmplitudePhotonBunch(self, diffraction_setup, incoming_bunch):
        """
        Calculates the diffraction/transmission given by the setup.
        :param diffraction_setup: The diffraction setup.
        :return: PhotonBunch object made up of diffracted/transmitted photons.
        """
        # Create PhotonBunch instance.
        outgoing_bunch = ComplexAmplitudePhotonBunch([])

        # Retrieve the photon bunch from the diffraction setup.
        # incoming_bunch = diffraction_setup.incomingPhotons()

        # Check that photon_bunch is indeed a PhotonBunch object.
        if not isinstance(incoming_bunch, ComplexAmplitudePhotonBunch):
            raise Exception("The incoming photon bunch must be a ComplexAmplitudePhotonBunch object!")

        # Raise calculation start.
        self._onCalculationStart()

        for index, polarized_photon in enumerate(incoming_bunch):

            # Raise OnProgress event if progressed by 10 percent.
            self._onProgressEveryTenPercent(index, len(incoming_bunch))

            outgoing_complex_amplitude_photon = self.calculateDiffractedComplexAmplitudePhoton(diffraction_setup,
                                                                        polarized_photon)
            # Add result of current deviation.
            outgoing_bunch.addPhoton(outgoing_complex_amplitude_photon)

        # Raise calculation end.
        self._onCalculationEnd()

        # Return diffraction results.
        return outgoing_bunch