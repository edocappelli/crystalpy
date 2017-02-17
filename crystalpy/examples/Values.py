"""
Allows the user to input the diffraction parameters or to use default values.
"""
from crystalpy.diffraction import GeometryType
from crystalpy.util.StokesVector import StokesVector


class Interval(object):

    def __init__(self, inf, sup):
        """
        Constructor.
        """
        self.inf = inf
        self.sup = sup

    def print(self):
        """
        Prints the interval to video.
        """
        print("\nInterval[{0}, {1}]\n".format(self.inf, self.sup))

    def check_extrema(self):
        """
        Returns True if inf <= sup. False if inf > sup.
        :return: bool
        """
        if self.inf > self.sup:
            return False

        return True

    def __eq__(self, candidate):
        """
        Returns True if self == candidate.
        """
        if (self.inf == candidate.inf and
                self.sup == candidate.sup):
            return True

        return False

    def __ne__(self, candidate):
        """
        Returns True if self != candidate.
        """
        return not self == candidate


class Values(object):

    def __init__(self):
        """
        Constructor.
        Prompt the user for values or use default ones.
        """
        while True:

            use_default = input("Do you want to use default values? (Y/N) ")

            if use_default in ("y", "Y", "Yes", "yes", "YES"):

                self._use_default()
                break

            elif use_default in ("n", "N", "No", "no", "NO"):

                use_hirano = input("Do you want to test one of the figures from Hirano et al.? (Y/N) ")

                if use_hirano in ("y", "Y", "Yes", "yes", "YES"):

                    self._use_hirano()
                    break

                else:
                    self._set_values()
                    break

            else:
                print("The input could not be interpreted. Try again.\n")

    def _use_default(self):
        """
        Use the default diffraction parameters.
        """
        # diffraction part
        self.geometry_type = GeometryType.BraggDiffraction()
        self.crystal_name = "Si"
        self.thickness = 0.01  # centimetres
        self.miller_h = 1  # int
        self.miller_k = 1  # int
        self.miller_l = 1  # int
        self.asymmetry_angle = 0.0  # degrees
        self.azimuthal_angle = 90.0  # degrees
        self.energy_min = 8.0  # keV
        self.energy_max = 8.0  # keV
        self.energy_points = 1  # int
        self.angle_deviation_min = -100  # micro radians
        self.angle_deviation_max = 100  # micro radians
        self.angle_deviation_points = 200  # int
        self.deg = True  # phase results in degrees by default.
        self.inclination_angle = 45.0  # degrees
        self.phase_inf_limit = -180  # degrees
        self.phase_sup_limit = 180  # degrees

        # polarization part
        self.incoming_stokes_vector = StokesVector([1, 1, 0, 0])

    def _set_values(self):
        """
        Set values for the diffraction parameters.
        """
        while True:

            self.geometry_type = int(input("\nBragg diffraction [0]\nBragg transmission [1]\n"
                                           "Laue diffraction [2]\nLaue transmission [3]\n\nChoose one type: "))
            if self.geometry_type == 0:
                self.geometry_type = GeometryType.BraggDiffraction()
                break
            elif self.geometry_type == 1:
                self.geometry_type = GeometryType.BraggTransmission()
                break
            elif self.geometry_type == 2:
                self.geometry_type = GeometryType.LaueDiffraction()
                break
            elif self.geometry_type == 3:
                self.geometry_type = GeometryType.LaueTransmission()
                break
            else:
                print("The input could not be interpreted. Try again.\n")

        while True:
            try:
                self.crystal_name = str(input("\nCrystal name: "))
                self.thickness = float(input("\nCrystal thickness [cm]: "))
                self.miller_h = int(input("\nMiller H: "))
                self.miller_k = int(input("\nMiller K: "))
                self.miller_l = int(input("\nMiller L: "))
                self.asymmetry_angle = float(input("\nAsymmetry angle [degrees]: "))
                self.azimuthal_angle = float(input("\nAzimuthal angle [degrees]: "))
                self.energy_min = float(input("\nMinimum energy [keV]: "))
                self.energy_max = float(input("\nMaximum energy [keV]: "))
                self.energy_points = int(input("\nNumber of energy points: "))
                self.angle_deviation_min = float(input("\nMinimum deviation from Bragg angle [micro radians]: "))
                self.angle_deviation_max = float(input("\nMaximum deviation from Bragg angle [micro radians]: "))
                self.angle_deviation_points = int(input("\nNumber of deviation points: "))
                self.inclination_angle = float(input("\nInclination angle [degrees]: "))
                self.deg = bool(int(input("\nShould the phase be represented in:\nradians[0]?\ndegrees[1]?")))
                self._stokes_parameters = list()
                self.phase_inf_limit = -180  # degrees
                self.phase_sup_limit = 180  # degrees

                for i in range(4):
                    new_element = float(input("\nStokes parameter S{}: ".format(i)))
                    self._stokes_parameters.append(new_element)
                self.incoming_stokes_vector = StokesVector(self._stokes_parameters)

                break

            except ValueError:
                print("\nAn error occurred. Try again.")

    def _use_hirano(self):
        """
        Use the settings for the figures in the Hirano et al. article.
        """
        print("\n----------\nUse the settings specified in the figures from\n"
              "K.Hirano et al., 'Perfect crystal X-ray phase retarders'\n----------\n")
        figure_number = int(input("Which figure do you want to reproduce? [e.g. for 'Fig.4' write '4'] "))

        # I call the files where I store the parameters "Fig_$_Hirano.dat", where $ = figure_number.
        file_name = "Fig_{}_Hirano.dat".format(figure_number)

        # Open file.
        file = open(file_name, "r")

        # Read file.
        file_data = list()
        for line in file:
            file_data.append(line.partition("  #")[0])

        self.geometry_type = str(file_data[0])

        if self.geometry_type == "BraggDiffraction":
            self.geometry_type = GeometryType.BraggDiffraction()

        elif self.geometry_type == "BraggTransmission":
            self.geometry_type = GeometryType.BraggTransmission()

        elif self.geometry_type == "LaueDiffraction":
            self.geometry_type = GeometryType.LaueDiffraction()

        elif self.geometry_type == "LaueTransmission":
            self.geometry_type = GeometryType.LaueTransmission()

        else:
            raise Exception("The file content couldn't be read.")

        try:
            self.crystal_name = str(file_data[1])
            self.thickness = float(file_data[2])
            self.miller_h = int(file_data[3])
            self.miller_k = int(file_data[4])
            self.miller_l = int(file_data[5])
            self.asymmetry_angle = float(file_data[6])
            self.energy_min = float(file_data[7])
            self.energy_max = float(file_data[8])
            self.energy_points = int(file_data[9])
            self.angle_deviation_min = float(file_data[10])
            self.angle_deviation_max = float(file_data[11])
            self.angle_deviation_points = int(file_data[12])
            self.inclination_angle = float(file_data[13])
            self.deg = bool(file_data[14])
            self._stokes_parameters = list()

            for i in range(4):
                new_element = float(file_data[15 + i])
                self._stokes_parameters.append(new_element)
            self.incoming_stokes_vector = StokesVector(self._stokes_parameters)

            # read the phase limitations.
            self.phase_inf_limit = float(file_data[19])
            self.phase_sup_limit = float(file_data[20])

            # read the intervals to set to zero.
            self.intervals_number = int(file_data[21])  # number of intervals.
            self.intervals = list()

            # if there are no intervals to set to zero close the file and return.
            if self.intervals_number == 0:
                file.close()
                return

            # else start reading the intervals.
            else:
                for i in range(self.intervals_number):

                    interval = Interval(float(file_data[22 + i * 2]), float(file_data[23 + i * 2]))
                    self.intervals.append(interval)

        except ValueError:
            raise Exception("The file content couldn't be read.")

        finally:
            file.close()

    def print(self):
        """
        Prints the values.
        """
        print("\nVALUES:"
              "\ngeometry type: {geometry_type}"
              "\ncrystal name: {crystal_name}"
              "\nthickness: {thickness}"
              "\nmiller H: {miller_h}"
              "\nmiller K: {miller_k}"
              "\nmiller L: {miller_l}"
              "\nasymmetry angle: {asymmetry_angle}"
              "\nazimuthal angle: {azimuthal_angle}"
              "\nenergy minimum: {energy_min}"
              "\nenergy maximum: {energy_max}"
              "\nenergy points: {energy_points}"
              "\ndeviation angle minimum: {angle_deviation_min}"
              "\ndeviation angle maximum: {angle_deviation_max}"
              "\ndeviation angle points: {angle_deviation_points}"
              "\ndegrees: {deg}"
              "\ninclination angle: {inclination_angle}"
              "\ninferior phase limit: {phase_inf_limit}"
              "\nsuperior phase limit: {phase_sup_limit}"
              "\nincoming Stokes vector: {incoming_stokes_vector}".format(geometry_type=self.geometry_type.description(),
                                                                          crystal_name=self.crystal_name,
                                                                          thickness=self.thickness,
                                                                          miller_h=self.miller_h,
                                                                          miller_k=self.miller_k,
                                                                          miller_l=self.miller_l,
                                                                          asymmetry_angle=self.asymmetry_angle,
                                                                          azimuthal_angle=self.azimuthal_angle,
                                                                          energy_min=self.energy_min,
                                                                          energy_max=self.energy_max,
                                                                          energy_points=self.energy_points,
                                                                          angle_deviation_min=self.angle_deviation_min,
                                                                          angle_deviation_max=self.angle_deviation_max,
                                                                          angle_deviation_points=self.angle_deviation_points,
                                                                          deg=self.deg,
                                                                          inclination_angle=self.inclination_angle,
                                                                          phase_inf_limit=self.phase_inf_limit,
                                                                          phase_sup_limit=self.phase_sup_limit,
                                                                          incoming_stokes_vector=self.incoming_stokes_vector.components()))
