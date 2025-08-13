import os
from pathlib import Path
from django.apps import apps
from django.test import TestCase

from epolyscat_django_app.trecx_utils import Linp

BASE_DIR = Path(__file__).parent


class LinpWithMasterFileTestCase(TestCase):
    def setUp(self) -> None:
        trecx_settings = apps.get_app_config("epolyscat_django_app").APPLICATION_SETTINGS[
            "tRecX"
        ]
        self.linp = Linp(trecx_settings["MASTER_LINP"])

    def test_absorption_axis(self):
        item = self.linp.item("Absorption", "axis", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "Rn")
        self.assertEqual(item.default, "NONE")
        self.assertEqual(item.docu, "absorption parameters for this axis")
        self.assertFalse(item.empty)

    def test_absorption_axis_empty(self):
        item = self.linp.item("Absorption", "axis", 2)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "EMPTY")
        self.assertEqual(item.default, "NONE")
        self.assertEqual(item.docu, "absorption parameters for this axis")
        self.assertTrue(item.empty)

    def test_absorption_axis_no_third_line(self):
        item = self.linp.item("Absorption", "axis", 3)
        self.assertIsNone(item)

    def test_absorption_kind(self):
        item = self.linp.item("Absorption", "kind", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "NOT_FOUND")
        self.assertEqual(item.default, "ECS")
        self.assertEqual(item.docu, "complex scaling kind")
        self.assertTrue(item.empty)

    def test_axis_alternate(self):
        item = self.linp.item("Axis", "alternate", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "NOT_FOUND")
        self.assertEqual(item.default, "")
        self.assertEqual(
            item.docu,
            'someName[order=10] will create same axis and FEs , but replace order by 10, store index as "someName"',
        )
        self.assertTrue(item.empty)

    def test_axis_functions_2(self):
        item = self.linp.item("Axis", "functions", 2)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "polExp[1.]")
        self.assertEqual(item.default, "polynomial")
        self.assertEqual(item.docu, "which functions to use")
        self.assertFalse(item.empty)

    def test_axis_functions_4(self):
        item = self.linp.item("Axis", "functions", 4)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "NOT_FOUND")
        self.assertEqual(item.default, "assocLegendre{Phi}")
        self.assertEqual(item.docu, "which functions to use")
        self.assertTrue(item.empty)

    def test_title(self):
        item = self.linp.item("Title", "", 1)
        self.assertIsNotNone(item)
        self.assertEqual(
            item.inputValue, "Tutorial example  - photo-electron spectra of hydrogen"
        )
        self.assertEqual(item.default, "tRecX calculation")
        self.assertEqual(item.docu, "title that will be printed to output")
        self.assertFalse(item.empty)

    def test_EXPERT_FDgrid(self):
        item = self.linp.item("_EXPERT_", "FDgrid", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "NOT_FOUND")
        self.assertEqual(item.default, "exp")
        self.assertEqual(item.docu, "kind of grid spacing: exp,sudden,smooth")
        self.assertTrue(item.empty)

    def test_Laser_lambda_at_I_2(self):
        item = self.linp.item("Laser", "lambda@I/2(nm)", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "NOT_FOUND")
        self.assertEqual(item.default, "20")
        self.assertEqual(
            item.docu, "linear chirp - wave length at half intensity after peak"
        )
        self.assertTrue(item.empty)

    def test_Laser_IWcm2(self):
        item = self.linp.item("Laser", "I(W/cm2)", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "1.e15")
        self.assertEqual(item.default, "NO_VALUE")
        self.assertEqual(item.docu, "peak intensity")
        self.assertFalse(item.empty)

    def test_Surface_definition(self):
        item = self.linp.item("Surface", "definition", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "NOT_FOUND")
        self.assertEqual(item.default, "not given")
        self.assertEqual(
            item.docu,
            "convert to: Spherical[R,Lmax{,Mmax}]...spherical basis with maximal angular momenta",
        )
        self.assertTrue(item.empty)

    def test_Spectrum_plot(self):
        item = self.linp.item("Spectrum", "plot", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "NOT_FOUND")
        self.assertEqual(item.default, "")
        self.assertEqual(
            item.docu,
            "(mutable) read from ampl-file and plot, default depends on coordinates:, total, partial, partialIon, intPhi...2d integrated over Phi, cutPhi[n]...n cuts at phi angles, cutEta[n]...n cuts at eta=cos(theta), sumZ...sum over Z (for XZ or YZ), sumY...sum over Y (for XY or YZ), sumX...sum over X (for XY or XZ), cone...integrate cone[phiC,etaC,gamC] axis phiC,etaC, opening gamC, zone...integrate zone[the0,the1] and phi in [0,2pi], partialRn1...sum over Rn2, partial waves, sumRn1...sum over Rn1, sumRn2...sum over Rn2, kZ...z-momentum, lPartial...partial waves summed over m, JAD[thMin{:thMax:n},eMin{:eMax:n}{~eV}{,phiMin{:phiMax:n}}]...joint angular distribution, user...read plot from Plot_user in input, grid[file]...draw the plot info from file, default (coordinate-dependent)",
        )
        self.assertTrue(item.empty)


class LinpWithTutorialFileTestCase(TestCase):
    def setUp(self) -> None:
        linp_file = BASE_DIR / "test-data" / "linp"
        with linp_file.open() as f:
            self.linp = Linp(f)

    def test_absorption_axis(self):
        item = self.linp.item("Absorption", "axis", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "NONE")
        self.assertEqual(item.default, "NONE")
        self.assertEqual(item.docu, "absorption parameters for this axis")
        # For the older style of linp that are without \inputValue{} tags,
        # there's no good way to determine if the item had no value, so empty
        # should always return true
        self.assertFalse(item.empty)

    def test_missing_absorption_axis_line(self):
        item = self.linp.item("Absorption", "axis", 2)
        self.assertIsNone(item)

    def test_absorption_kind(self):
        item = self.linp.item("Absorption", "kind", 1)
        self.assertIsNone(item)

    def test_axis_upper_end(self):
        item = self.linp.item("Axis", "upper end", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "10")
        self.assertEqual(item.default, "1.7976931348623e+308")
        self.assertEqual(
            item.docu,
            "upper boundary of the axis",
        )
        self.assertFalse(item.empty)

    def test_axis_functions(self):
        item = self.linp.item("Axis", "functions", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "polynomial")
        self.assertEqual(item.default, "polynomial")
        self.assertEqual(item.docu, "which functions to use")
        self.assertFalse(item.empty)

    def test_axis_functions_4(self):
        item = self.linp.item("Axis", "functions", 4)
        self.assertIsNone(item)

    def test_title(self):
        item = self.linp.item("Title", "", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "Tutorial example  - 1d harmonic oscillator")
        self.assertEqual(item.default, "tRecX calculation")
        self.assertEqual(item.docu, "title that will be printed to output")
        self.assertFalse(item.empty)

    def test_EXPERT_FDgrid(self):
        item = self.linp.item("_EXPERT_", "FDgrid", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "exp")
        self.assertEqual(item.default, "exp")
        self.assertEqual(item.docu, "kind of grid spacing: exp,sudden,smooth")
        self.assertFalse(item.empty)

    def test_Operator_hamiltonian(self):
        item = self.linp.item("Operator", "hamiltonian", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "0.5<d_1_d>+0.5<Q*Q>")
        self.assertEqual(item.default, "NO_VALUE")
        self.assertEqual(item.docu, "Hamiltonian definition")
        self.assertFalse(item.empty)

    def test_Flag_generate_linp(self):
        item = self.linp.item("Flag", "generate-linp", 0)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "FLAG_ONLY")
        self.assertEqual(item.default, "FLAG_ONLY")
        self.assertEqual(
            item.docu,
            "only generate list-of-inputs file",
        )
        self.assertFalse(item.empty)

    def test_Spectrum_plot(self):
        item = self.linp.item("Spectrum", "plot", 1)
        self.assertIsNotNone(item)
        self.assertEqual(item.inputValue, "")
        self.assertEqual(item.default, "")
        self.assertEqual(
            item.docu,
            "(mutable) read from ampl-file and plot, default depends on coordinates:, total, partial, partialIon, intPhi...2d integrated over Phi, cutPhi[n]...n cuts at phi angles, cutEta[n]...n cuts at eta=cos(theta), sumZ...sum over Z (for XZ or YZ), sumY...sum over Y (for XY or YZ), sumX...sum over X (for XY or XZ), cone...integrate cone[phiC,etaC,gamC] axis phiC,etaC, opening gamC, zone...integrate zone[the0,the1] and phi in [0,2pi], partialRn1...sum over Rn2, partial waves, sumRn1...sum over Rn1, sumRn2...sum over Rn2, kZ...z-momentum, lPartial...partial waves summed over m, JAD[thMin{:thMax:n},eMin{:eMax:n}{~eV}{,phiMin{:phiMax:n}}]...joint angular distribution, user...read plot from Plot_user in input, grid[file]...draw the plot info from file, default (coordinate-dependent)",
        )
        self.assertFalse(item.empty)
