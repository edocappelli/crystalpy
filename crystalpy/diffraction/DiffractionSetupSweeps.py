"""
Represents a diffraction setup (from where it inherits) and includes a photon beam
with scanning ephoton energy and deviation angle (scattering plane is YZ plane)

Units are in SI except for photon energy in eV. Angles in radians.

"""
import numpy

from crystalpy.diffraction.DiffractionSetup import DiffractionSetup

from crystalpy.util.Photon import Photon
from crystalpy.util.PhotonBunch import PhotonBunch


class DiffractionSetupSweeps(DiffractionSetup):

    def __init__(self,
                 geometry_type, crystal_name, thickness,
                 miller_h, miller_k, miller_l,
                 asymmetry_angle,
                 azimuthal_angle,
                 energy_min,
                 energy_max,
                 energy_points,
                 angle_deviation_min,
                 angle_deviation_max,
                 angle_deviation_points):
        """
        Constructor.
        :param geometry_type: GeometryType (BraggDiffraction,...).
        :param crystal_name: The name of the crystal, e.g. Si.
        :param thickness: The crystal thickness.
        :param miller_h: Miller index H.
        :param miller_k: Miller index K.
        :param miller_l: Miller index L.
        :param asymmetry_angle: The asymmetry angle between surface normal and Bragg normal.
        :param azimuthal_angle: The angle between the projection of the Bragg normal
                                on the crystal surface plane and the x axis.
        :param energy_min: The minimum energy.
        :param energy_max: The maximum energy.
        :param energy_points: Number of energy points.
        :param angle_deviation_min: Minimal angle deviation.
        :param angle_deviation_max: Maximal angle deviation.
        :param angle_deviation_points: Number of deviations points.
        """
        energies = numpy.linspace(energy_min,
                               energy_max,
                               energy_points)

        deviations = numpy.linspace(angle_deviation_min,
                                 angle_deviation_max,
                                 angle_deviation_points)

        # Create an "info setup" solely for the determination of deviation angles.
        info_setup = DiffractionSetup(geometry_type=geometry_type,
                                      crystal_name=crystal_name,
                                      thickness=thickness,
                                      miller_h=miller_h,
                                      miller_k=miller_k,
                                      miller_l=miller_l,
                                      asymmetry_angle=asymmetry_angle,
                                      azimuthal_angle=azimuthal_angle,)
                                      # incoming_photons=[Photon(energy_min, Vector(0, 0, -1))])


        # Call base constructor.
        DiffractionSetup.__init__(self,
                                  geometry_type=geometry_type,
                                  crystal_name=crystal_name,
                                  thickness=thickness,
                                  miller_h=miller_h,
                                  miller_k=miller_k,
                                  miller_l=miller_l,
                                  asymmetry_angle=asymmetry_angle,
                                  azimuthal_angle=azimuthal_angle,)

                                  # incoming_photons=photons)


        # Create photons according to sweeps.
        photons = PhotonBunch() # list()
        for energy in energies:
            for deviation in deviations:
                direction = info_setup.incomingPhotonDirection(energy, deviation)
                incoming_photon = Photon(energy, direction)
                # photons.append(incoming_photon)
                photons.addPhoton(incoming_photon)

        self._incoming_photons = photons
        # srio@esrf.eu: in theory, not needed as this info is in _incoming_photons
        # but buffering this accelerates a lot the calculations
        self._deviations = None
        self._energies = None


    #
    # overwritten methods
    #

    def toDictionary(self):
        """
        Returns this setup in InfoDictionary form.
        :return: InfoDictionary form of this setup.
        """
        info_dict = DiffractionSetup.toDictionary(self)

        info_dict["Minimum energy"] = str(self.energyMin())
        info_dict["Maximum energy"] = str(self.energyMax())
        info_dict["Number of energy points"] = str(self.energyPoints())
        info_dict["Angle deviation minimum"] = "%.2e" % (self.angleDeviationMin())
        info_dict["Angle deviation maximum"] = "%.2e" % (self.angleDeviationMax())
        info_dict["Angle deviation points"] = str(self.angleDeviationPoints())

        return info_dict

    def __eq__(self, candidate):
        """
        Determines if two setups are equal.
        :param candidate: Instance to compare to.
        :return: True if the two instances are equal. False otherwise.
        """
        if self._geometry_type != candidate.geometryType():
            return False

        if self._crystal_name != candidate.crystalName():
            return False

        if self._thickness != candidate.thickness():
            return False

        if self._miller_h != candidate.millerH():
            return False

        if self._miller_k != candidate.millerK():
            return False

        if self._miller_l != candidate.millerL():
            return False

        if self._asymmetry_angle != candidate.asymmetryAngle():
            return False

        if self._azimuthal_angle != candidate.azimuthalAngle():
            return False

        if self.energyMin() != candidate.energyMin():
            return False

        if self.energyMax() != candidate.energyMax():
            return False

        if self.energyPoints() != candidate.energyPoints():
            return False

        if self.angleDeviationMin() != candidate.angleDeviationMin():
            return False

        if self.angleDeviationMax() != candidate.angleDeviationMax():
            return False

        if self.angleDeviationPoints() != candidate.angleDeviationPoints():
            return False

        # All members are equal so are the instances.
        return True

    #
    # end overwritten methods
    #

    def incomingPhotons(self):
        """
        Returns the incoming photons.
        :return: A list of photons.
        """
        #TODO return the PhotonBunch object ?
        return self._incoming_photons.getListOfPhotons()

    def energyMin(self):
        """
        Returns the minimum energy in eV.
        :return: The minimum energy in eV.
        """
        return self.energies().min()

    def energyMax(self):
        """
        Returns the maximum energy in eV.
        :return: The maximum energy in eV.
        """
        return self.energies().max()

    def energyPoints(self):
        """
        Returns the number of energy points.
        :return: Number of energy points.
        """
        return self.energies().shape[0]

    def energies(self):
        """
        Returns the energies of this setup.
        :return: The angle deviations grid.
        """
        if self._energies is None:
            self._energies = numpy.unique(numpy.array([photon.energy() for photon in self._incoming_photons.getListOfPhotons()]))

        return self._energies

    def angleDeviationMin(self):
        """
        Returns the minimal angle deviation.
        :return: Minimal angle deviation.
        """
        return self.angleDeviationGrid().min()

    def angleDeviationMax(self):
        """
        Returns the maximal angle deviation.
        :return: Maximal angle deviation.
        """
        return self.angleDeviationGrid().max()

    def angleDeviationPoints(self):
        """
        Returns the angle deviation points.
        :return: Angle deviation points.
        """
        return self.angleDeviationGrid().shape[0]

    def angleDeviationGrid(self):
        """
        Returns the grid of angle deviations according to this setup.
        :return: The angle deviations grid.
        """
        if self._deviations is None:
            self._deviations = numpy.array([self.deviationOfIncomingPhoton(photon) for photon in self._incoming_photons.getListOfPhotons()])


        return self._deviations

