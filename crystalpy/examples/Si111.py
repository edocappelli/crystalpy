import numpy as np

from crystalpy.diffraction.GeometryType import GeometryType, BraggDiffraction
from crystalpy.diffraction.DiffractionSetupSweeps import DiffractionSetupSweeps
from crystalpy.diffraction.Diffraction import Diffraction

from crystalpy.PlotData1D import PlotData1D
from crystalpy.polarization.MuellerDiffraction import MuellerDiffraction

import matplotlib.pyplot as plt
from crystalpy.Values import Values


from crystalpy.polarization.StokesVector import StokesVector



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
    info_dict = result.diffractionSetup().asInfoDictionary()
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
    info_dict = result.diffraction_setup.asInfoDictionary()
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

# def four_plots(x1,y1,x2,y2,x3,y3,x4,y4,title="",xtitle="",ytitle="",xrange=None,yrange=None,show=True):
#     """
#     Creates four plots in a window
#
#     :param x1: abscissas for plot 1
#     :param y1: ordinates for plot 1
#     :param x2: abscissas for plot 2
#     :param y2: ordinates for plot 2
#     :param x3: abscissas for plot 3
#     :param y3: ordinates for plot 3
#     :param x4: abscissas for plot 4
#     :param y4: ordinates for plot 4
#     :param title: a string or list of 4 strings with title
#     :param xtitle: a string or list of 4 strings with title for X
#     :param ytitle: a string or list of 4 strings with title for Y
#     :param xrange: the X range for all plots
#     :param yrange: the Y range for all plots
#     :param show:
#     :return:
#     """
#
#     if isinstance(title,list):
#         Title = title
#     else:
#         Title = [title,title,title,title]
#
#     if isinstance(xtitle,list):
#         Xtitle = xtitle
#     else:
#         Xtitle = [xtitle,xtitle,xtitle,xtitle]
#
#     if isinstance(ytitle,list):
#         Ytitle = ytitle
#     else:
#         Ytitle = [ytitle,ytitle,ytitle,ytitle]
#
#     # Create subplots.
#     f, ((ax00, ax01), (ax10, ax11)) = plt.subplots(2, 2, sharex="all", sharey="all")
#
#     ax00.plot(x1,y1, "-")
#     ax00.set_title(  Title[0])
#     ax00.set_xlabel(Xtitle[0])
#     ax00.set_ylabel(Ytitle[0])
#     ax00.set_xlim(xrange)
#     ax00.set_ylim(yrange)
#
#
#     ax01.plot(x2,y2, "-")
#     ax01.set_title(  Title[1])
#     ax01.set_xlabel(Xtitle[1])
#     ax01.set_ylabel(Ytitle[1])
#     ax01.set_xlim(xrange)
#     ax01.set_ylim(yrange)
#
#
#     ax10.plot(x3,y3, "-")
#     ax10.set_title(  Title[2])
#     ax10.set_xlabel(Xtitle[2])
#     ax10.set_ylabel(Ytitle[2])
#     ax10.set_xlim(xrange)
#     ax10.set_ylim(yrange)
#
#
#     ax11.plot(x4,y4, "-")
#     ax11.set_title(  Title[3])
#     ax11.set_xlabel(Xtitle[3])
#     ax11.set_ylabel(Ytitle[3])
#     ax11.set_xlim(xrange)
#     ax11.set_ylim(yrange)
#
#     if show: plt.show()

if __name__ == "__main__":

    # Get values from user or use default values.
    # values = Values()
    # values.print()

    # Create a diffraction setup.
    # At this stage I translate angles in radians, energy in eV and all other values in SI units.
    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetupSweeps(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 1e-3,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0.0,                              # radians
                                               azimuthal_angle        = 0.0,                              # radians
                                               energy_min             = 8000.0,                           # eV
                                               energy_max             = 8000.0,                           # eV
                                               energy_points          = 1,                                # int
                                               angle_deviation_min    = -100e-6,                          # radians
                                               angle_deviation_max    = 100e-6,                           # radians
                                               angle_deviation_points = 500)                              # int

    # Create a Diffraction object.
    diffraction = Diffraction()

    # Create a DiffractionResult object holding the results of the diffraction calculations.
    print("\nCalculating the diffraction results...")
    diffraction_result = diffraction.calculateDiffraction(diffraction_setup)


    #
    # plots
    #
    photon_energies = diffraction_result.energies()
    deviation_angles = diffraction_result.angleDeviations()
    print("Number of energy points: %d"%photon_energies.size)
    print("Number of angular points: %d"%deviation_angles.size)

    print("_intensity shape: ",diffraction_result._intensities.shape)
    print("_phases shape: ",diffraction_result._phases.shape)

    from srxraylib.plot.gol import plot, four_plots
    # plot( 1e6*deviation_angles,diffraction_result._intensities[0,:,0],
    #       1e6*deviation_angles,diffraction_result._intensities[0,:,1],
    #       1e6*deviation_angles,diffraction_result._intensities[0,:,2],
    #       title="Intensity for photon energy = %4.3f "%photon_energies[0],
    #       xtitle="Deviation angle urad",ytitle="Reflectivity",
    #       legend=['s-pol','p-pol','p/s ratio',],show=False)
    #
    # plot( 1e6*deviation_angles,diffraction_result._phases[0,:,0],
    #       1e6*deviation_angles,diffraction_result._phases[0,:,1],
    #       1e6*deviation_angles,diffraction_result._phases[0,:,2],
    #       title="Phase for photon energy = %4.3f "%photon_energies[0],
    #       xtitle="Deviation angle urad",ytitle="Reflectivity",
    #       legend=['s-pol','p-pol','p minus s pol'],show=False)


    # Create a PlotData1D object.
    print("\nCreating the diffraction profile plots...")

    # values = Values()
    # plot_1d = plot_diffraction_1d(diffraction_result, values.deg)

    # if False:
    #     # Unwrap the phases.
    #     print("\nUnwrapping the phase data...")
    #     phase_limits = (values.phase_inf_limit, values.phase_sup_limit)
    #     plot_1d[3].smart_unwrap(values.intervals, values.intervals_number, phase_limits, values.deg)
    #     plot_1d[4].smart_unwrap(values.intervals, values.intervals_number, phase_limits, values.deg)
    #     plot_1d[5].smart_unwrap(values.intervals, values.intervals_number, phase_limits, values.deg)

    # # Plot the diffraction results.
    # intensity_phase_plot(plot_1d, values)
    # #
    # # # Create a plot following the representation used  in:
    # # # K.Hirano et al., 'Perfect Crystal X-ray phase retarders' (1993).
    # hirano_plot(plot_1d, values)


    # Create a MuellerDiffraction object.
    mueller_diffraction = MuellerDiffraction(diffraction_result,
                                             StokesVector([1,1,0,1]),
                                             inclination_angle=0) # np.pi*45/180)

    # Create a MullerResult object.
    print("\nCalculating the Stokes vector...")
    mueller_result = mueller_diffraction.calculate_stokes()

    # # Create a PlotData1D object.
    # print("\nCreating the Stokes parameters plots...")
    # plot_1d = plot_stokes_1d(mueller_result)

    # Plot the Stokes vectors.
    # stokes_plot(plot_1d, values)

    four_plots(1e6*deviation_angles,mueller_result._s0[0],
               1e6*deviation_angles,mueller_result._s1[0],
               1e6*deviation_angles,mueller_result._s2[0],
               1e6*deviation_angles,mueller_result._s3[0],
               title=["S0","S1","S2","S3"],xtitle="Deviation angle [urad]",
               yrange=[-1,1])

    # Plot the degree of circular polarization.
    # polarization_degree_plot(plot_1d, values)

    plt.show()
