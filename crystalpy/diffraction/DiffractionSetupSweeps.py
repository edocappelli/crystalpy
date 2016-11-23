"""
Represents a diffraction setup.
Except for energy all units are in SI. Energy is in eV. Angles in radians.
"""
import numpy as np

from orangecontrib.crystal.diffraction.DiffractionSetup import DiffractionSetup

from orangecontrib.crystal.util.Photon import Photon
from orangecontrib.crystal.util.Vector import Vector


class DiffractionSetupSweeps(DiffractionSetup):

    def __init__(self, geometry_type, crystal_name, thickness,
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
        energies = np.linspace(energy_min,
                               energy_max,
                               energy_points)

        deviations = np.linspace(angle_deviation_min,
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
                                      azimuthal_angle=azimuthal_angle,
                                      incoming_photons=[Photon(energy_min, Vector(0, 0, -1))])

        # Create photons according to sweeps.
        photons = list()
        for energy in energies:
            for deviation in deviations:
                direction = info_setup.incomingPhotonDirection(energy, deviation)
                incoming_photon = Photon(energy, direction)
                photons.append(incoming_photon)

        # Call base constructor.
        DiffractionSetup.__init__(self,
                                  geometry_type=geometry_type,
                                  crystal_name=crystal_name,
                                  thickness=thickness,
                                  miller_h=miller_h,
                                  miller_k=miller_k,
                                  miller_l=miller_l,
                                  asymmetry_angle=asymmetry_angle,
                                  azimuthal_angle=azimuthal_angle,
                                  incoming_photons=photons)
