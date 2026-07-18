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


def test_analysis_selection_preserves_order_and_rejects_duplicates():
    script = textwrap.dedent(
        """
        const fs = require("fs");
        const path = require("path");
        const babel = require("@babel/core");

        const sourcePath = path.join(
          process.cwd(), "src", "utils", "workflow-analysis-selection.js"
        );
        const source = fs.readFileSync(sourcePath, "utf8");
        const compiled = babel.transformSync(source, {
          plugins: ["@babel/plugin-transform-modules-commonjs"],
        }).code;
        const moduleObject = { exports: {} };
        new Function("module", "exports", compiled)(moduleObject, moduleObject.exports);

        const {
          addAnalysisApplication,
          normalizeAnalysisApplications,
        } = moduleObject.exports;
        const allowed = ["CnvMath", "CnvMatLab", "CnvLinFull"];
        const normalized = normalizeAnalysisApplications(
          ["CnvLinFull", "CnvMath", "CnvLinFull", "Unknown"],
          allowed,
          "CnvMath",
        );
        const added = addAnalysisApplication(normalized, "CnvMatLab", allowed);
        const duplicate = addAnalysisApplication(added, "CnvMath", allowed);

        if (JSON.stringify(normalized) !== JSON.stringify(["CnvLinFull", "CnvMath"])) {
          throw new Error(`Unexpected normalized order: ${JSON.stringify(normalized)}`);
        }
        if (JSON.stringify(added) !== JSON.stringify(["CnvLinFull", "CnvMath", "CnvMatLab"])) {
          throw new Error(`Unexpected added order: ${JSON.stringify(added)}`);
        }
        if (JSON.stringify(duplicate) !== JSON.stringify(added)) {
          throw new Error("Duplicate utility was added");
        }
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_analysis_selection_reorders_and_keeps_one_required_utility():
    script = textwrap.dedent(
        """
        const fs = require("fs");
        const path = require("path");
        const babel = require("@babel/core");

        const sourcePath = path.join(
          process.cwd(), "src", "utils", "workflow-analysis-selection.js"
        );
        const source = fs.readFileSync(sourcePath, "utf8");
        const compiled = babel.transformSync(source, {
          plugins: ["@babel/plugin-transform-modules-commonjs"],
        }).code;
        const moduleObject = { exports: {} };
        new Function("module", "exports", compiled)(moduleObject, moduleObject.exports);

        const {
          moveAnalysisApplication,
          removeAnalysisApplication,
        } = moduleObject.exports;
        const original = ["CnvMath", "CnvMatLab", "CnvLinFull"];
        const moved = moveAnalysisApplication(original, 2, -1);
        const boundary = moveAnalysisApplication(moved, 0, -1);
        const removed = removeAnalysisApplication(moved, "CnvMatLab");
        const retained = removeAnalysisApplication(["CnvMath"], "CnvMath");

        if (JSON.stringify(moved) !== JSON.stringify(["CnvMath", "CnvLinFull", "CnvMatLab"])) {
          throw new Error(`Unexpected move result: ${JSON.stringify(moved)}`);
        }
        if (JSON.stringify(boundary) !== JSON.stringify(moved)) {
          throw new Error("Boundary move changed the list");
        }
        if (JSON.stringify(removed) !== JSON.stringify(["CnvMath", "CnvLinFull"])) {
          throw new Error(`Unexpected remove result: ${JSON.stringify(removed)}`);
        }
        if (JSON.stringify(retained) !== JSON.stringify(["CnvMath"])) {
          throw new Error("The final required utility was removed");
        }
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr
