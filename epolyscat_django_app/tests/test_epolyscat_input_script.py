import subprocess
import textwrap

from pathlib import Path


APP_DIR = Path(__file__).resolve().parents[1]


def _run_node_script(script):
    return subprocess.run(
        ["node", "-e", script],
        cwd=APP_DIR,
        text=True,
        capture_output=True,
        check=False,
    )


def test_table_updates_patch_uploaded_input_file_without_dropping_unknown_records():
    script = textwrap.dedent(
        """
        const fs = require("fs");
        const path = require("path");
        const babel = require("@babel/core");

        const sourcePath = path.join(process.cwd(), "src", "utils", "epolyscat-input-script.js");
        const source = fs.readFileSync(sourcePath, "utf8");
        const compiled = babel.transformSync(source, {
          plugins: ["@babel/plugin-transform-modules-commonjs"],
        }).code;
        const moduleObject = { exports: {} };
        new Function("module", "exports", compiled)(moduleObject, moduleObject.exports);

        const { patchEPolyScatInputScript } = moduleObject.exports;
        const original = [
          "# original title",
          "LMax 4",
          "CustomLegacyRecord keep me",
          "FileName 'PlotData' 'old-cross-sections.dat' 'REWIND'",
        ].join("\\n");
        const outputDefinitions = [
          {
            fileType: "PlotData",
            valueKey: "plotDataFile",
            disposition: "REWIND",
          },
        ];
        const patched = patchEPolyScatInputScript(original, {
          title: "patched title",
          lMax: 9,
          eMax: 20,
          plotDataFile: "new-cross-sections.dat",
        }, outputDefinitions);

        function assertContains(value) {
          if (!patched.includes(value)) {
            throw new Error(`Expected patched script to contain: ${value}\\n${patched}`);
          }
        }

        function assertNotContains(value) {
          if (patched.includes(value)) {
            throw new Error(`Expected patched script not to contain: ${value}\\n${patched}`);
          }
        }

        assertContains("# patched title");
        assertContains("LMax 9");
        assertContains("EMax 20");
        assertContains("CustomLegacyRecord keep me");
        assertContains("FileName 'PlotData' 'new-cross-sections.dat' 'REWIND'");
        assertNotContains("old-cross-sections.dat");
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_table_field_patch_does_not_append_unmodified_default_records():
    script = textwrap.dedent(
        """
        const fs = require("fs");
        const path = require("path");
        const babel = require("@babel/core");

        const sourcePath = path.join(process.cwd(), "src", "utils", "epolyscat-input-script.js");
        const source = fs.readFileSync(sourcePath, "utf8");
        const compiled = babel.transformSync(source, {
          plugins: ["@babel/plugin-transform-modules-commonjs"],
        }).code;
        const moduleObject = { exports: {} };
        new Function("module", "exports", compiled)(moduleObject, moduleObject.exports);

        const { patchEPolyScatInputScript } = moduleObject.exports;
        const original = [
          "# original title",
          "LMax 4",
          "CustomLegacyRecord keep me",
        ].join("\\n");
        const patched = patchEPolyScatInputScript(original, {
          lMax: 11,
          eMax: 20,
        }, [], { changedKeys: ["lMax"] });

        if (!patched.includes("LMax 11")) {
          throw new Error(`Expected LMax to be patched\\n${patched}`);
        }
        if (!patched.includes("CustomLegacyRecord keep me")) {
          throw new Error(`Expected unknown record to be preserved\\n${patched}`);
        }
        if (patched.includes("EMax 20")) {
          throw new Error(`Expected unmodified EMax not to be appended\\n${patched}`);
        }
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_parse_uploaded_test03_input_preserves_engform_and_fegeeng_records():
    script = textwrap.dedent(
        """
        const fs = require("fs");
        const path = require("path");
        const babel = require("@babel/core");

        const sourcePath = path.join(process.cwd(), "src", "utils", "epolyscat-input-script.js");
        const source = fs.readFileSync(sourcePath, "utf8");
        const compiled = babel.transformSync(source, {
          plugins: ["@babel/plugin-transform-modules-commonjs"],
        }).code;
        const moduleObject = { exports: {} };
        new Function("module", "exports", compiled)(moduleObject, moduleObject.exports);

        const { parseEPolyScatInputScript } = moduleObject.exports;
        const contents = [
          "# input file for test03",
          "# electron scattering from N2 molden SCF, DCS calculation",
          "  LMax   15     # maximum l to be used for wave functions",
          "  EMax  50.0    # EMax, maximum asymptotic energy in eV",
          "  EngForm      # Energy formulas",
          "   0 0         # charge, formula type",
          "  FegeEng 13.0   # Energy correction (in eV) used in the fege potential",
          "  ScatContSym 'SG'  # Scattering symmetry",
          "  LMaxK    4     # Maximum l in the K matirx",
          "  ScatEng 3.0 4.0 5.0 6.0",
          "Convert '$pt/test03.molden2012' 'molden'",
          "GetBlms",
          "ExpOrb",
          "GetPot",
          "GrnType 1",
          "  ScatContSym 'B2U' # Scattering symmetry",
          "Scat",
          "TotalCrossSection",
          "EDCS",
        ].join("\\n");

        const parsed = parseEPolyScatInputScript(contents);

        function assertEqual(name, actual, expected) {
          if (actual !== expected) {
            throw new Error(`${name}: expected ${expected}, got ${actual}\\n${JSON.stringify(parsed, null, 2)}`);
          }
        }

        assertEqual("lMax", parsed.lMax, 15);
        assertEqual("eMax", parsed.eMax, 50);
        assertEqual("engFormCharge", parsed.engFormCharge, 0);
        assertEqual("engFormType", parsed.engFormType, 0);
        if ("engFormTerms" in parsed) {
          throw new Error(`Expected missing EngForm term count not to be parsed\\n${JSON.stringify(parsed, null, 2)}`);
        }
        assertEqual("fegeEng", parsed.fegeEng, 13);
        assertEqual("lMaxK", parsed.lMaxK, 4);
        assertEqual("scatEng", parsed.scatEng, "3.0 4.0 5.0 6.0");
        assertEqual("convertSource", parsed.convertSource, "$pt/test03.molden2012");
        assertEqual("convertFormat", parsed.convertFormat, "molden");
        assertEqual("scatContSym", parsed.scatContSym, "B2U");
        assertEqual("calculationKind", parsed.calculationKind, "scattering");
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_parse_uploaded_input_omits_fields_not_present_in_file():
    script = textwrap.dedent(
        """
        const fs = require("fs");
        const path = require("path");
        const babel = require("@babel/core");

        const sourcePath = path.join(process.cwd(), "src", "utils", "epolyscat-input-script.js");
        const source = fs.readFileSync(sourcePath, "utf8");
        const compiled = babel.transformSync(source, {
          plugins: ["@babel/plugin-transform-modules-commonjs"],
        }).code;
        const moduleObject = { exports: {} };
        new Function("module", "exports", compiled)(moduleObject, moduleObject.exports);

        const { parseEPolyScatInputScript } = moduleObject.exports;
        const parsed = parseEPolyScatInputScript([
          "LMax 15",
          "EngForm",
          " 0 0",
          "FegeEng 13.0",
        ].join("\\n"));

        if (parsed.lMax !== 15 || parsed.engFormCharge !== 0 || parsed.engFormType !== 0 || parsed.fegeEng !== 13) {
          throw new Error(`Expected present records to parse\\n${JSON.stringify(parsed, null, 2)}`);
        }

        ["lMaxI", "eMax", "engFormTerms", "scatEng", "convertSource"].forEach(key => {
          if (key in parsed) {
            throw new Error(`Expected ${key} to be omitted when absent\\n${JSON.stringify(parsed, null, 2)}`);
          }
        });
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_rtf_input_contents_are_normalized_before_parsing():
    script = textwrap.dedent(
        """
        const fs = require("fs");
        const path = require("path");
        const babel = require("@babel/core");

        const sourcePath = path.join(process.cwd(), "src", "utils", "epolyscat-input-script.js");
        const source = fs.readFileSync(sourcePath, "utf8");
        const compiled = babel.transformSync(source, {
          plugins: ["@babel/plugin-transform-modules-commonjs"],
        }).code;
        const moduleObject = { exports: {} };
        new Function("module", "exports", compiled)(moduleObject, moduleObject.exports);

        const {
          normalizeEPolyScatInputContents,
          parseEPolyScatInputScript,
        } = moduleObject.exports;

        const rtf = [
          "{\\\\rtf1\\\\ansi\\\\ansicpg936{\\\\fonttbl\\\\f0\\\\fmodern Courier;}",
          "\\\\pard\\\\pardeftab720",
          "\\\\f0\\\\fs26 \\\\cf0 \\\\strokec2 #\\\\ input file for test03\\\\par",
          "LMax 15\\\\par",
          "EMax 50.0\\\\par",
          "FegeEng 13.0}",
        ].join("\\n");
        const normalized = normalizeEPolyScatInputContents(rtf);

        if (normalized.includes("\\\\rtf") || normalized.includes("fonttbl") || normalized.includes("strokec2")) {
          throw new Error(`Expected RTF controls to be removed\\n${normalized}`);
        }
        if (!normalized.includes("# input file for test03") || !normalized.includes("LMax 15")) {
          throw new Error(`Expected plain ePolyScat text\\n${normalized}`);
        }

        const parsed = parseEPolyScatInputScript(normalized);
        if (parsed.lMax !== 15 || parsed.eMax !== 50 || parsed.fegeEng !== 13) {
          throw new Error(`Expected normalized RTF to parse\\n${JSON.stringify(parsed, null, 2)}\\n${normalized}`);
        }
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_parser_removes_rtf_escape_backslashes_from_title_and_scateng():
    script = textwrap.dedent(
        """
        const fs = require("fs");
        const path = require("path");
        const babel = require("@babel/core");

        const sourcePath = path.join(process.cwd(), "src", "utils", "epolyscat-input-script.js");
        const source = fs.readFileSync(sourcePath, "utf8");
        const compiled = babel.transformSync(source, {
          plugins: ["@babel/plugin-transform-modules-commonjs"],
        }).code;
        const moduleObject = { exports: {} };
        new Function("module", "exports", compiled)(moduleObject, moduleObject.exports);

        const {
          normalizeEPolyScatInputContents,
          parseEPolyScatInputScript,
        } = moduleObject.exports;

        const contents = [
          "#\\\\ electron scattering from N2 molden SCF, DCS calculation",
          "ScatEng 3.0 4.0 5.0 6.0\\\\",
        ].join("\\n");
        const normalized = normalizeEPolyScatInputContents(contents);
        const parsed = parseEPolyScatInputScript(normalized);

        if (parsed.title !== "electron scattering from N2 molden SCF, DCS calculation") {
          throw new Error(`Expected title without backslash, got ${parsed.title}\\n${normalized}`);
        }
        if (parsed.scatEng !== "3.0 4.0 5.0 6.0") {
          throw new Error(`Expected ScatEng without trailing backslash, got ${parsed.scatEng}\\n${normalized}`);
        }
        if (normalized.includes("\\\\")) {
          throw new Error(`Expected normalized content not to retain RTF escape slashes\\n${normalized}`);
        }
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_workflow_output_binding_uses_manual_specific_file_flow():
    script = textwrap.dedent(
        """
        const fs = require("fs");
        const path = require("path");
        const babel = require("@babel/core");

        const sourcePath = path.join(process.cwd(), "src", "utils", "workflow-file-linking.js");
        const source = fs.readFileSync(sourcePath, "utf8");
        const compiled = babel.transformSync(source, {
          plugins: ["@babel/plugin-transform-modules-commonjs"],
        }).code;
        const moduleObject = { exports: {} };
        new Function("module", "exports", compiled)(moduleObject, moduleObject.exports);

        const { buildWorkflowOutputInputBinding } = moduleObject.exports;

        function assertEqual(name, actual, expected) {
          if (actual !== expected) {
            throw new Error(`${name}: expected ${expected}, got ${actual}`);
          }
        }

        const dataGenerationOutputs = [
          { name: "molcas.log", dataProductURI: "airavata-dp://molcas-log" },
          { name: "molden.dat", dataProductURI: "airavata-dp://molden" },
          { name: "gaussian.log", dataProductURI: "airavata-dp://gaussian-log" },
        ];

        const gaussianBinding = buildWorkflowOutputInputBinding({
          outputFiles: dataGenerationOutputs,
          targetStageId: "ePolyScat_Run",
          sourceApplicationId: "Gaussian16",
        });
        assertEqual("Gaussian input", gaussianBinding.outputFile.name, "gaussian.log");
        assertEqual("Gaussian convert source", gaussianBinding.dataEntryValues.convertSource, "gaussian.log");
        assertEqual("Gaussian convert format", gaussianBinding.dataEntryValues.convertFormat, "gaussian");

        const molcasBinding = buildWorkflowOutputInputBinding({
          outputFiles: dataGenerationOutputs,
          targetStageId: "ePolyScat_Run",
          sourceApplicationId: "OpenMolcas",
        });
        assertEqual("OpenMolcas input", molcasBinding.outputFile.name, "molden.dat");
        assertEqual("OpenMolcas convert source", molcasBinding.dataEntryValues.convertSource, "molden.dat");
        assertEqual("OpenMolcas convert format", molcasBinding.dataEntryValues.convertFormat, "molden");

        const noCnvLinFullBinding = buildWorkflowOutputInputBinding({
          outputFiles: [
            { name: "cross-sections.dat", dataProductURI: "airavata-dp://plot" },
            { name: "matrix-elements.idy", dataProductURI: "airavata-dp://idy" },
          ],
          targetStageId: "Analysis",
          targetApplicationId: "CnvLinFull",
          requiredFileName: "DumpOut",
        });
        if (noCnvLinFullBinding !== null) {
          throw new Error(`Expected CnvLinFull not to bind a generic dat file`);
        }

        const cnvLinFullBinding = buildWorkflowOutputInputBinding({
          outputFiles: [
            { name: "matrix-elements.idy", dataProductURI: "airavata-dp://idy" },
            { name: "test15dumpidy.dat", dataProductURI: "airavata-dp://dumpidy" },
          ],
          targetStageId: "Analysis",
          targetApplicationId: "CnvLinFull",
          requiredFileName: "DumpOut",
        });
        assertEqual("CnvLinFull input", cnvLinFullBinding.outputFile.name, "test15dumpidy.dat");
        assertEqual("CnvLinFull input name", cnvLinFullBinding.inputFileName, "DumpOut");
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr
