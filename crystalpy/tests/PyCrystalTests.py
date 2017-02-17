import unittest

# diffraction
from crystalpy.tests.diffraction.GeometryTypeTest import GeometryTypeTest
from crystalpy.tests.diffraction.ComplexAmplitudeTest import ComplexAmplitudeTest
from crystalpy.tests.diffraction.PerfectCrystalDiffractionTest import PerfectCrystalDiffractionTest
from crystalpy.tests.diffraction.DiffractionSetupTest import DiffractionSetupTest
from crystalpy.tests.diffraction.DiffractionSetupSweepsTest import DiffractionSetupSweepsTest
from crystalpy.tests.diffraction.DiffractionTest import DiffractionTest
from crystalpy.tests.diffraction.DiffractionResultTest import DiffractionResultTest

# polarization
from crystalpy.tests.polarization.CrystalPhasePlateTest import CrystalPhasePlateTest
from crystalpy.tests.polarization.MuellerDiffractionTest import MuellerDiffractionTest
from crystalpy.tests.polarization.MuellerMatrixTest import MuellerMatrixTest
from crystalpy.tests.polarization.MuellerResultTest import MuellerResultTest
from crystalpy.tests.util.StokesVectorTest import StokesVectorTest

# util
from crystalpy.tests.util.VectorTest import VectorTest
from crystalpy.tests.util.PhotonTest import PhotonTest
from crystalpy.tests.util.PhotonBunchTest import PhotonBunchTest

from crystalpy.tests.util.PolarizedPhotonTest import PolarizedPhotonTest
from crystalpy.tests.util.PolarizedPhotonBunchTest import PolarizedPhotonBunchTest


def suite():
    """
    Gathers all the tests in a test suite.
    """
    suites = (
        # tests by Mark Glass.
        unittest.makeSuite(VectorTest, 'test'),
        unittest.makeSuite(PhotonTest, 'test'),
        unittest.makeSuite(PhotonBunchTest, 'test'),

        unittest.makeSuite(ComplexAmplitudeTest, 'test'),
        unittest.makeSuite(GeometryTypeTest, 'test'),
        unittest.makeSuite(PerfectCrystalDiffractionTest, 'test'),
        unittest.makeSuite(DiffractionSetupTest, 'test'),
        unittest.makeSuite(DiffractionSetupSweepsTest, 'test'),
        unittest.makeSuite(DiffractionTest, 'test'),
        unittest.makeSuite(DiffractionResultTest, 'test'),

        # tests by Edoardo.

        unittest.makeSuite(StokesVectorTest, "test"),
        unittest.makeSuite(PolarizedPhotonTest, "test"),
        unittest.makeSuite(PolarizedPhotonBunchTest, "test")   ,
        unittest.makeSuite(CrystalPhasePlateTest, "test"),
        unittest.makeSuite(MuellerDiffractionTest, "test"),
        unittest.makeSuite(MuellerMatrixTest, "test"),
        unittest.makeSuite(MuellerResultTest, "test"),

    )
    return unittest.TestSuite(suites)


if __name__ == "__main__":

    unittest.main(defaultTest="suite")
