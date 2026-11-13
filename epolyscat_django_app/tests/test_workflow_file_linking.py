import subprocess
import textwrap

from pathlib import Path

import pytest


APP_DIR = Path(__file__).resolve().parents[1]


def _run_node_script(script):
    return subprocess.run(
        ["node", "-e", script],
        cwd=APP_DIR,
        text=True,
        capture_output=True,
        check=False,
    )


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

        const {
          buildEPolyScatInputDataSelectionValues,
          buildWorkflowOutputInputBinding,
        } = moduleObject.exports;

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
        assertEqual("Gaussian convert source", gaussianBinding.dataEntryValues.convertSource, "$pt/gaussian.log");
        assertEqual("Gaussian convert format", gaussianBinding.dataEntryValues.convertFormat, "gaussian");

        const molcasBinding = buildWorkflowOutputInputBinding({
          outputFiles: dataGenerationOutputs,
          targetStageId: "ePolyScat_Run",
          sourceApplicationId: "OpenMolcas",
        });
        assertEqual("OpenMolcas input", molcasBinding.outputFile.name, "molden.dat");
        assertEqual("OpenMolcas convert source", molcasBinding.dataEntryValues.convertSource, "$pt/molden.dat");
        assertEqual("OpenMolcas convert format", molcasBinding.dataEntryValues.convertFormat, "molden");

        const uploadedMoldenValues = buildEPolyScatInputDataSelectionValues({
          name: "/Users/researcher/inputs/ch4-target.molden",
        });
        assertEqual("Uploaded Molden source", uploadedMoldenValues.convertSource, "$pt/ch4-target.molden");
        assertEqual("Uploaded Molden format", uploadedMoldenValues.convertFormat, "molden");

        const uploadedGaussianValues = buildEPolyScatInputDataSelectionValues({
          name: "methane-gaussian.log",
        });
        assertEqual("Uploaded Gaussian source", uploadedGaussianValues.convertSource, "$pt/methane-gaussian.log");
        assertEqual("Uploaded Gaussian format", uploadedGaussianValues.convertFormat, "gaussian");

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
