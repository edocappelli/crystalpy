"""
---OK---
"""
from collections import OrderedDict
import copy

import numpy as np

from crystalpy.examples.Values import Interval


class PlotData1D(object):
    """
    Represents a 1D plot. The graph data together with related information.
    """
    def __init__(self, title, title_x_axis, title_y_axis):
        """
        Constructor.
        :param title: Plot title.
        :param title_x_axis: X axis' title.
        :param title_y_axis: Y axis' title.
        """
        # Set titles.
        self.title = title
        self.title_x_axis = title_x_axis
        self.title_y_axis = title_y_axis

        # Initialize X and Y ranges.
        self.x_min = None
        self.x_max = None
        self.y_min = None
        self.y_max = None

        # Initialize X and Y data.
        self.x = None
        self.y = None

        # Initialize plot information to empty ordered dictionary.
        self._plot_info = OrderedDict()

    def set_x_min(self, x_min):
        """
        Sets x range minimum.
        :param x_min: X range minimum.
        """
        self.x_min = x_min

    def set_x_max(self, x_max):
        """
        Sets X range maximum.
        :param x_max: X range maximum.
        """
        self.x_max = x_max

    def set_y_min(self, y_min):
        """
        Sets Y range minimum.
        :param y_min: Y range minimum.
        """
        self.y_min = y_min

    def set_y_max(self, y_max):
        """
        Sets Y range maximum.
        :param y_max: Y range maximum.
        """
        self.y_max = y_max

    def set_x(self, x):
        """
        Sets X data.
        :param x: x data.
        """
        self.x = x

    def set_y(self, y):
        """
        Sets Y data.
        :param y: y data.
        """
        self.y = y

    def _set_interval_to_zero(self, indices, lower=True, upper=True):
        """
        Sets the y's to zero in certain intervals of x's (extrema included).
        :param indices: pair with the two extrema of the x interval.
        :param lower: if True include the lower end of the interval.
        :param upper: if True include the upper end of the interval.
        """
        try:
            inf_index = indices.inf
            sup_index = indices.sup

            # adjust the indices according to the lower and upper parameters.
            if not lower:
                inf_index += 1

            if not upper:
                sup_index -= 1

            # in the index range defined by inf_index and sup_index, set the y's to zero.
            for i in range(inf_index, sup_index + 1):
                self.y[i] = 0

        except TypeError:
            print("\nERROR: could not set the values to zero in the specified intervals.\n")

    def _unwrap_interval(self, indices, deg, lower=True, upper=True):
        """
        Unwraps the y data vector in a certain interval.
        :param indices: indices determining the interval to unwrap.
        :param deg: True if values are in degrees. False if radians.
        :param lower: if True include the lower end of the interval.
        :param upper: if True include the upper end of the interval.
        """
        inf_index = indices.inf
        sup_index = indices.sup

        # adjust the indices according to the lower and upper parameters.
        if not lower:
            inf_index += 1

        if not upper:
            sup_index -= 1

        # numpy.unwrap works on data in radians, so if the data is in degrees, it needs to be converted.
        if deg:
            self.y = np.deg2rad(self.y)

            # cut out the part to unwrap and then stitch it back on.
            temp = self.y[inf_index:sup_index + 1]
            self.y[inf_index:sup_index + 1] = np.unwrap(temp)

            # convert back to degrees.
            self.y = np.rad2deg(self.y)
            return

        # cut out the part to unwrap and then stitch it back on.
        temp = self.y[inf_index:sup_index + 1]
        self.y[inf_index:sup_index + 1] = np.unwrap(temp)

    def _optimize_interval(self, indices, phase_limits):
        """
        Takes an interval and restricts it so that the extrema match the points where the phase
        becomes bigger(smaller) than some upper(lower) limit.
        :param indices: indices corresponding to the interval to be optimized.
        :param phase_limits: the limits of the phase to be used for the optimization, [min, max].
        :return: indices of the optimized interval.
        """
        inf = indices.inf
        sup = indices.sup

        # check the intervals.
        if (self.y[inf] > phase_limits[1] or
                self.y[inf] < phase_limits[0]):
            print("\nERROR in PlotData1D._optimize_interval: First value in the interval exceeds limitations.")
            return indices

        if (self.y[sup] > phase_limits[1] or
                self.y[sup] < phase_limits[0]):
            print("\nERROR in PlotData1D._optimize_interval: Last value in the interval exceeds limitations.")
            return indices

        # starting from the lower end.
        i = inf  # counter initialization.
        while phase_limits[0] < self.y[i] < phase_limits[1]:
            i += 1

        # if the conditions are not satisfied for index i:
        new_inf = i - 1

        # starting from the upper end.
        i = sup  # counter initialization.
        while phase_limits[0] < self.y[i] < phase_limits[1]:
            i -= 1

        # if the conditions are not satisfied for index i:
        new_sup = i + 1

        new_indices = Interval(new_inf, new_sup)

        # check that the inf is smaller than (or equal to) the sup.
        if not new_indices.check_extrema():
            print("\nERROR in PlotData1D._optimize_interval: The phase might be undersampled.")
            return indices

        return new_indices

    def smart_unwrap(self, intervals, intervals_number, phase_limits, deg):
        """
        Unwraps data correctly by avoiding discontinuities.
        :param intervals: list of pairs. Each element is a pair with the two extrema of the x interval.
        :param phase_limits: min and max tolerable values for the phase plot, [min, max].
        :param intervals_number: number of intervals to set to zero.
        :param deg: True if values are in degrees. False if radians.
        """
        if intervals_number == 0:
            if deg:
                self.y = np.deg2rad(self.y)  # unwrap works with radians.
                self.y = np.unwrap(self.y)
                self.y = np.rad2deg(self.y)  # convert back to degrees.
                return

            self.y = np.unwrap(self.y)
            return

        # transform self.x into a numpy.ndarray object.
        x = np.asarray(self.x)

        # careful! only works with monotonic sequences.
        temp_index = x.argmin()

        for interval in intervals:
            inf = interval.inf
            sup = interval.sup

            # find the indices of the y array corresponding to inf and sup.
            inf_index = abs(x - inf).argmin()
            sup_index = abs(x - sup).argmin()

            # optimize the interval.
            indices = Interval(inf_index, sup_index)
            new_indices = self._optimize_interval(indices, phase_limits)

            # unwrap the data before the interval.
            indices_to_unwrap = Interval(temp_index, new_indices.inf)
            self._unwrap_interval(indices_to_unwrap, deg, lower=True, upper=False)

            # set the interval to zero.
            indices_to_set = new_indices
            self._set_interval_to_zero(indices_to_set, lower=True, upper=False)

            temp_index = new_indices.sup

        # careful! only works with monotonic sequences.
        indices_to_unwrap = Interval(temp_index, x.argmax())
        self._unwrap_interval(indices_to_unwrap, deg, lower=True, upper=True)

    def add_xy_point(self, x_point, y_point):
        """
        Adds an x-y point.
        :param x_point: x coordinate.
        :param y_point: y coordinate.
        """
        self.x.append(x_point)
        self.y.append(y_point)
        
    def add_plot_info(self, name, info):
        """
        Adds a plot info.
        :param name: Name of the info.
        :param info: The info.
        """
        self._plot_info[name] = info
        
    def plot_info(self):
        """
        Returns the plot info copy.
        :return: The plot info.
        """
        return copy.deepcopy(self._plot_info)
