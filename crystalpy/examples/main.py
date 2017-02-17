import numpy as np
import matplotlib.pyplot as plt

from crystalpy.diffraction.DiffractionSetupSweeps import DiffractionSetupSweeps
from crystalpy.diffraction.Diffraction import Diffraction
from crystalpy.examples.Values import Values
from crystalpy.examples.PlotData1D import PlotData1D
from crystalpy.polarization.MuellerDiffraction import MuellerDiffraction


def intensity_phase_plot(plot_1d, values):
    """
    Plot the diffraction results.
    :param plot_1d: PlotData1D object.
    :param values: Values object.
    """
    # Create subplots.
    f, ((ax1_intensity, ax2_intensity, ax3_intensity),
        (ax1_phase, ax2_phase, ax3_phase)) = plt.subplots(2, 3, sharex="all", sharey="row")

    ax1_intensity.plot(plot_1d[0].x, plot_1d[0].y, "b-")
    ax1_intensity.set_title(plot_1d[0].title)
    ax1_intensity.set_xlabel(plot_1d[0].title_x_axis)
    ax1_intensity.set_ylabel(plot_1d[0].title_y_axis)
    ax1_intensity.set_xlim([values.angle_deviation_min, values.angle_deviation_max])
    ax1_intensity.set_ylim([0, 1])

    ax1_phase.plot(plot_1d[3].x, plot_1d[3].y, "g-.")
    ax1_phase.set_title(plot_1d[3].title)
    ax1_phase.set_xlabel(plot_1d[3].title_x_axis)
    ax1_phase.set_ylabel(plot_1d[3].title_y_axis)
    ax1_phase.set_xlim([values.angle_deviation_min, values.angle_deviation_max])
    ax1_phase.set_ylim([values.phase_inf_limit, values.phase_sup_limit])

    ax2_intensity.plot(plot_1d[1].x, plot_1d[1].y, "b-")
    ax2_intensity.set_title(plot_1d[1].title)
    ax2_intensity.set_xlabel(plot_1d[1].title_x_axis)
    ax2_intensity.set_ylabel(plot_1d[1].title_y_axis)
    ax2_intensity.set_xlim([values.angle_deviation_min, values.angle_deviation_max])
    ax2_intensity.set_ylim([0, 1])

    ax2_phase.plot(plot_1d[4].x, plot_1d[4].y, "g-.")
    ax2_phase.set_title(plot_1d[4].title)
    ax2_phase.set_xlabel(plot_1d[4].title_x_axis)
    ax2_phase.set_ylabel(plot_1d[4].title_y_axis)
    ax2_phase.set_xlim([values.angle_deviation_min, values.angle_deviation_max])
    ax2_phase.set_ylim([values.phase_inf_limit, values.phase_sup_limit])

    ax3_intensity.plot(plot_1d[2].x, plot_1d[2].y, "b-")
    ax3_intensity.set_title(plot_1d[2].title)
    ax3_intensity.set_xlabel(plot_1d[2].title_x_axis)
    ax3_intensity.set_ylabel(plot_1d[2].title_y_axis)
    ax3_intensity.set_xlim([values.angle_deviation_min, values.angle_deviation_max])
    ax3_intensity.set_ylim([0, 1])

    ax3_phase.plot(plot_1d[5].x, plot_1d[5].y, "g-.")
    ax3_phase.set_title(plot_1d[5].title)
    ax3_phase.set_xlabel(plot_1d[5].title_x_axis)
    ax3_phase.set_ylabel(plot_1d[5].title_y_axis)
    ax3_phase.set_xlim([values.angle_deviation_min, values.angle_deviation_max])
    ax3_phase.set_ylim([values.phase_inf_limit, values.phase_sup_limit])


def hirano_plot(plot_1d, values):
    """
    Create a plot following the representation used  in:
    K.Hirano et al., 'Perfect Crystal X-ray phase retarders' (1993).
    :param plot_1d: PlotData1D object.
    :param values: Values object.
    """
    # Create subplots.
    f, ax = plt.subplots(1, 1)

    # Intensity plots.
    ax.plot(plot_1d[0].x, plot_1d[0].y, "b-", label="sigma intensity")
    ax.plot(plot_1d[1].x, plot_1d[1].y, "k--", label="pi intensity")
    plt.legend(loc="upper left")
    ax.set_title("Hirano plot")
    ax.set_xlabel(plot_1d[0].title_x_axis)
    ax.set_ylabel(plot_1d[0].title_y_axis)
    ax.set_xlim([values.angle_deviation_min, values.angle_deviation_max])

    # Phase difference plot.
    ax_bis = ax.twinx()  # phase and intensities share the same x axis.
    ax_bis.set_xlim([values.angle_deviation_min, values.angle_deviation_max])
    ax_bis.set_ylim([values.phase_inf_limit, values.phase_sup_limit])
    ax_bis.plot(plot_1d[5].x, plot_1d[5].y, "g-.", label="phase retardation")
    ax_bis.set_ylabel(plot_1d[5].title_y_axis)
    plt.legend(loc="center left")


def stokes_plot(plot_1d, values):
    """
    Plot the Stokes vectors.
    :param plot_1d: PlotData1D object.
    :param values: Values object.
    """
    # Create subplots.
    f, ((ax00, ax01), (ax10, ax11)) = plt.subplots(2, 2, sharex="all", sharey="all")

    ax00.plot(plot_1d[0].x, plot_1d[0].y, "-")
    ax00.set_title(plot_1d[0].title)
    ax00.set_xlabel(plot_1d[0].title_x_axis)
    ax00.set_ylabel(plot_1d[0].title_y_axis)
    ax00.set_xlim([values.angle_deviation_min, values.angle_deviation_max])
    ax00.set_ylim([-1.0, 1.0])

    ax01.plot(plot_1d[1].x, plot_1d[1].y, "-")
    ax01.set_title(plot_1d[1].title)
    ax01.set_xlabel(plot_1d[1].title_x_axis)
    ax01.set_ylabel(plot_1d[1].title_y_axis)
    ax01.set_xlim([values.angle_deviation_min, values.angle_deviation_max])
    ax01.set_ylim([-1.0, 1.0])

    ax10.plot(plot_1d[2].x, plot_1d[2].y, "-")
    ax10.set_title(plot_1d[2].title)
    ax10.set_xlabel(plot_1d[2].title_x_axis)
    ax10.set_ylabel(plot_1d[2].title_y_axis)
    ax10.set_xlim([values.angle_deviation_min, values.angle_deviation_max])
    ax10.set_ylim([-1.0, 1.0])

    ax11.plot(plot_1d[3].x, plot_1d[3].y, "-")
    ax11.set_title(plot_1d[3].title)
    ax11.set_xlabel(plot_1d[3].title_x_axis)
    ax11.set_ylabel(plot_1d[3].title_y_axis)
    ax11.set_xlim([values.angle_deviation_min, values.angle_deviation_max])
    ax11.set_ylim([-1.0, 1.0])


def polarization_degree_plot(plot_1d, values):
    """
    Plot the degree of circular polarization.
    :param plot_1d: PlotData1D object.
    :param values: Values object.
    """
    f, polarization_degree = plt.subplots(1, 1)
    polarization_degree.plot(plot_1d[4].x, plot_1d[4].y, "b-")
    polarization_degree.set_title(plot_1d[4].title)
    polarization_degree.set_xlabel(plot_1d[4].title_x_axis)
    polarization_degree.set_ylabel(plot_1d[4].title_y_axis)
    polarization_degree.set_xlim([values.angle_deviation_min, values.angle_deviation_max])
    polarization_degree.set_ylim([-1, 1])


def plot_diffraction_1d(result, deg):
    """
    Returns this result instance in PlotData1D representation.
    :param deg: if False the phase is expressed in radians, if True in degrees.
    """
    # Distinguish between the strings "phase in deg" and "phase in rad".
    if deg:
        phase_string = "Phase in deg"
    else:
        phase_string = "Phase in rad"

    # Retrieve setup information.
    info_dict = result.diffractionSetup().toDictionary()
    info_dict["Bragg angle"] = str(result.braggAngle())

    # Retrieve angles of the results.
    angles_in_um = [i * 1e+6 for i in result.angleDeviations()]

    # Define inner function to duplicate info for every plot.
    def addPlotInfo(info_dict, energy, angles_in_um, data):
        plot_data = PlotData1D(data[0], data[1], data[2])
        plot_data.set_x(angles_in_um)
        plot_data.set_y(data[3])
        for key, value in info_dict.items():
            plot_data.add_plot_info(key, value)
        plot_data.add_plot_info("Energy", str(energy))
        return plot_data

    plots = []
    for energy in result.energies():
        # Intensity S polarization.
        categories = []

        s_intensity = ("Intensity - Polarization S",
                       "Angle deviation in urad",
                       "Intensity",
                       result.sIntensityByEnergy(energy))
        plots.append(addPlotInfo(info_dict, energy, angles_in_um, s_intensity))

        p_intensity = ("Intensity - Polarization P",
                       "Angle deviation in urad",
                       "Intensity",
                       result.pIntensityByEnergy(energy))
        plots.append(addPlotInfo(info_dict, energy, angles_in_um, p_intensity))

        intensity_difference = ("Intensity difference",
                                "Angle deviation in urad",
                                "Intensity",
                                result.differenceIntensityByEnergy(energy))
        plots.append(addPlotInfo(info_dict, energy, angles_in_um, intensity_difference))

        s_phase = ("Phase - Polarization S",
                   "Angle deviation in urad",
                   phase_string,
                   result.sPhaseByEnergy(energy, deg))
        plots.append(addPlotInfo(info_dict, energy, angles_in_um, s_phase))

        p_phase = ("Phase - Polarization P",
                   "Angle deviation in urad",
                   phase_string,
                   result.pPhaseByEnergy(energy, deg))
        plots.append(addPlotInfo(info_dict, energy, angles_in_um, p_phase))

        phase_difference = ("Phase difference",
                            "Angle deviation in urad",
                            phase_string,
                            result.differencePhaseByEnergy(energy, deg))
        plots.append(addPlotInfo(info_dict, energy, angles_in_um, phase_difference))

    return plots


def plot_stokes_1d(result):
    """
    Returns this result instance in PlotData1D representation.
    """
    # Retrieve setup information.
    info_dict = result.diffraction_setup.toDictionary()
    info_dict["Bragg angle"] = str(result.diffraction_result.braggAngle())

    # Retrieve angles of the results.
    angles_in_urad = [i * 1e+6 for i in result.angle_deviations()]

    # Define inner function to duplicate info for every plot.
    def add_all_plot_info(info_dict, energy, angles_in_urad, data):
        plot_data = PlotData1D(data[0], data[1], data[2])
        plot_data.set_x(angles_in_urad)
        plot_data.set_y(data[3])
        for key, value in info_dict.items():  # dict.items() returns a list of (key, value) tuples
            plot_data.add_plot_info(key, value)
        plot_data.add_plot_info("Energy", str(energy))
        return plot_data

    plots = list()
    for energy in result.energies():
        s0 = ("Stokes parameter S0",
              "Angle deviation in urad",
              "S0",
              result.s0_by_energy(energy))
        plots.append(add_all_plot_info(info_dict, energy, angles_in_urad, s0))

        s1 = ("Stokes parameter S1",
              "Angle deviation in urad",
              "S1",
              result.s1_by_energy(energy))
        plots.append(add_all_plot_info(info_dict, energy, angles_in_urad, s1))

        s2 = ("Stokes parameter S2",
              "Angle deviation in urad",
              "S2",
              result.s2_by_energy(energy))
        plots.append(add_all_plot_info(info_dict, energy, angles_in_urad, s2))

        s3 = ("Stokes parameter S3",
              "Angle deviation in urad",
              "S3",
              result.s3_by_energy(energy))
        plots.append(add_all_plot_info(info_dict, energy, angles_in_urad, s3))

        polarization_degree = ("Degree of circular polarization",
                               "Angle deviation in urad",
                               "Circular polarization",
                               result.polarization_degree_by_energy(energy))
        plots.append(add_all_plot_info(info_dict, energy, angles_in_urad, polarization_degree))

    return plots

if __name__ == "__main__":

    # Get values from user or use default values.
    values = Values()
    values.print()

    # Create a diffraction setup.
    # At this stage I translate angles in radians, energy in eV and all other values in SI units.
    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetupSweeps(geometry_type=values.geometry_type,         # GeometryType object
                                               crystal_name=values.crystal_name,           # string
                                               thickness=float(values.thickness) * 1e-2,   # meters
                                               miller_h=int(values.miller_h),              # int
                                               miller_k=int(values.miller_k),              # int
                                               miller_l=int(values.miller_l),              # int
                                               asymmetry_angle=float(values.asymmetry_angle) / 180 * np.pi,   # radians
                                               azimuthal_angle=float(values.azimuthal_angle) / 180 * np.pi,   # radians
                                               energy_min=float(values.energy_min) * 1e3,  # eV
                                               energy_max=float(values.energy_max) * 1e3,  # eV
                                               energy_points=int(values.energy_points),    # int
                                               angle_deviation_min=float(values.angle_deviation_min) * 1e-6,  # radians
                                               angle_deviation_max=float(values.angle_deviation_max) * 1e-6,  # radians
                                               angle_deviation_points=int(values.angle_deviation_points))     # int

    # Create a Diffraction object.
    diffraction = Diffraction()

    # Create a DiffractionResult object holding the results of the diffraction calculations.
    print("\nCalculating the diffraction results...")
    diffraction_result = diffraction.calculateDiffraction(diffraction_setup)

    # Create a PlotData1D object.
    print("\nCreating the diffraction profile plots...")
    plot_1d = plot_diffraction_1d(diffraction_result, values.deg)

    if False:
        # Unwrap the phases.
        print("\nUnwrapping the phase data...")
        phase_limits = (values.phase_inf_limit, values.phase_sup_limit)
        plot_1d[3].smart_unwrap(values.intervals, values.intervals_number, phase_limits, values.deg)
        plot_1d[4].smart_unwrap(values.intervals, values.intervals_number, phase_limits, values.deg)
        plot_1d[5].smart_unwrap(values.intervals, values.intervals_number, phase_limits, values.deg)

    # Plot the diffraction results.
    intensity_phase_plot(plot_1d, values)

    # Create a plot following the representation used  in:
    # K.Hirano et al., 'Perfect Crystal X-ray phase retarders' (1993).
    hirano_plot(plot_1d, values)

    # Create a MuellerDiffraction object.
    mueller_diffraction = MuellerDiffraction(diffraction_result,
                                             values.incoming_stokes_vector,
                                             values.inclination_angle*np.pi/180.0)

    # Create a MullerResult object.
    print("\nCalculating the Stokes vector...")
    mueller_result = mueller_diffraction.calculate_stokes()

    # Create a PlotData1D object.
    print("\nCreating the Stokes parameters plots...")
    plot_1d = plot_stokes_1d(mueller_result)

    # Plot the Stokes vectors.
    stokes_plot(plot_1d, values)

    # Plot the degree of circular polarization.
    polarization_degree_plot(plot_1d, values)

    plt.show()
