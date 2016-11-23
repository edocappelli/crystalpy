"""
Exception classes.
"""


class DiffractionException(Exception):
    def __init__(self, exception_text):
        super(DiffractionException, self).__init__(exception_text)


class ReflectionImpossibleException(DiffractionException):
    def __init__(self):
        super(ReflectionImpossibleException, self).__init__("Impossible geometry. "
                                                            "Asymmetry angle larger than Bragg angle in Bragg geometry. "
                                                            "No reflection possible.")


class TransmissionImpossibleException(DiffractionException):
    def __init__(self):
        super(TransmissionImpossibleException, self).__init__("Impossible geometry. "
                                                              "Asymmetry angle smaller than Bragg angle in Laue geometry. "
                                                              "No transmission possible.")


class StructureFactorF0isZeroException(DiffractionException):
    def __init__(self):
        super(StructureFactorF0isZeroException, self).__init__("Structure factor for F_0 is zero.")


class StructureFactorFHisZeroException(DiffractionException):
    def __init__(self):
        super(StructureFactorFHisZeroException, self).__init__("Structure factor for H=(hkl) is zero. "
                                                               "Forbidden reflection for given Miller indices?")


class StructureFactorFHbarIsZeroException(DiffractionException):
    def __init__(self):
        super(StructureFactorFHbarIsZeroException, self).__init__("Structure factor for H_bar=(-h,-k,-l) is zero. "
                                                                  "Forbidden reflection for given Miller indices?")
