from crystalpy.diffraction.GeometryType import BraggDiffraction
from crystalpy.diffraction.DiffractionSetupSweeps import DiffractionSetupSweeps
from crystalpy.diffraction.Diffraction import Diffraction

from crystalpy.polarization.MuellerDiffraction import MuellerDiffraction
from crystalpy.util.StokesVector import StokesVector


def calculate_standard_interface():

    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetupSweeps(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 1e-2



                                               ,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
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
    # Now the Mueller/Stokes calculation from the diffraction results
    #

    mueller_diffraction = MuellerDiffraction(diffraction_result,
                                             StokesVector([1,0,1,0]),
                                             inclination_angle=0.0) #np.pi*45/180)

    # Create a MullerResult object.
    print("\nCalculating the Stokes vector...")
    mueller_result = mueller_diffraction.calculate_stokes()

    return mueller_result


def make_plots(mueller_result):
    #
    # plots
    #
    diffraction_result = mueller_result.diffraction_result

    photon_energies = diffraction_result.energies()
    deviation_angles = diffraction_result.angleDeviations()

    print("Number of energy points: %d"%photon_energies.size)
    print("Number of angular points: %d"%deviation_angles.size)
    print("_intensity shape: ",diffraction_result._intensities.shape)
    print("_phases shape: ",diffraction_result._phases.shape)

    from srxraylib.plot.gol import plot, four_plots
    plot( 1e6*deviation_angles,diffraction_result._intensities[0,:,0],
          1e6*deviation_angles,diffraction_result._intensities[0,:,1],
          1e6*deviation_angles,diffraction_result._intensities[0,:,2],
          title="Intensity for photon energy = %4.3f "%photon_energies[0],
          xtitle="Deviation angle urad",ytitle="Reflectivity",
          legend=['s-pol','p-pol','p/s ratio',],show=False)

    plot( 1e6*deviation_angles,diffraction_result._phases[0,:,0],
          1e6*deviation_angles,diffraction_result._phases[0,:,1],
          1e6*deviation_angles,diffraction_result._phases[0,:,2],
          title="Phase for photon energy = %4.3f "%photon_energies[0],
          xtitle="Deviation angle urad",ytitle="Reflectivity",
          legend=['s-pol','p-pol','p minus s pol'],show=False)

    # Stokes
    four_plots(1e6*deviation_angles,mueller_result._s0[0],
               1e6*deviation_angles,mueller_result._s1[0],
               1e6*deviation_angles,mueller_result._s2[0],
               1e6*deviation_angles,mueller_result._s3[0],
               title=["S0","S1","S2","S3"],xtitle="Deviation angle [urad]",
               yrange=[-1,1],show=False)

    # Plot the degree of circular polarization.
    plot(1e6*deviation_angles,mueller_result._s3[0]/mueller_result._s0[0],yrange=[-1,1],
         title="Circular Polarization S3/S0",xtitle="Deviation angle [urad]",ytitle="S3/S0",show=True)




if __name__ == "__main__":
    mueller_result = calculate_standard_interface()
    make_plots(mueller_result)