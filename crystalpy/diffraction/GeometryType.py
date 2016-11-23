"""
Represents geometry types/setups: Bragg diffraction, BraggTransmission, Laue diffraction, Laue transmission.
"""


class GeometryType(object):
    def __init__(self, description):
        """
        Constructor.
        :param description: Description of the geometry type, e.g. "Bragg transmission"
        """
        self._description = description

    def description(self):
        """
        Returns the description of this geometry type.
        :return: Description of this geometry type.
        """
        return self._description

    def __eq__(self, candidate):
        """
        Determines if two instances are equal.
        :param candidate: Instances to compare to.
        :return: True if both instances are equal. Otherwise False.
        """
        return self.description() == candidate.description()

    def __ne__(self, candidate):
        """
        Determines if two instances are not equal.
        :param candidate: Instances to compare.
        :return: True if both instances are not equal. Otherwise False.
        """
        return not self == candidate

    def __hash__(self):
        """
        Returns the hash value of this instance.
        :return: Hash value of this instance.
        """
        # As hash value just use the hash of the description.
        return hash(self._description)

    @staticmethod
    def allGeometryTypes():
        """
        Returns all possible geometry types.
        :return: All possible geometry types.
        """
        return [BraggDiffraction(),
                LaueDiffraction(),
                BraggTransmission(),
                LaueTransmission()]


class LaueDiffraction(GeometryType):
    """
    Represents Laue diffraction.
    """
    def __init__(self):
        super(LaueDiffraction, self).__init__("Laue diffraction")


class BraggDiffraction(GeometryType):
    """
    Represents Bragg diffraction.
    """
    def __init__(self):
        super(BraggDiffraction, self).__init__("Bragg diffraction")


class LaueTransmission(GeometryType):
    """
    Represents Laue transmission.
    """
    def __init__(self):
        super(LaueTransmission, self).__init__("Laue transmission")


class BraggTransmission(GeometryType):
    """
    Represents Bragg transmission.
    """
    def __init__(self):
        super(BraggTransmission, self).__init__("Bragg transmission")
