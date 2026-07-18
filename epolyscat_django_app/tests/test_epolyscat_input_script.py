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


def test_lossless_document_round_trip_preserves_order_comments_unknown_records_and_line_endings():
    script = textwrap.dedent(
        r'''
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
          parseEPolyScatDocument,
          serializeEPolyScatDocument,
        } = moduleObject.exports;
        const original = [
          "# repeated records are executable state\r\n",
          "ScatContSym 'SG'  # first calculation\r\n",
          "Scat\r\n",
          "CustomFutureRecord keep this exactly\r\n",
          "EngForm\r\n",
          " 0 2\r\n",
          " 2\r\n",
          " 2.0 -1.0 0\r\n",
          " 1.0  0.5 1\r\n",
          "ScatContSym 'SU'\r\n",
          "Scat\r\n",
        ].join("");

        const document = parseEPolyScatDocument(original);
        const serialized = serializeEPolyScatDocument(document);
        if (serialized !== original) {
          throw new Error(`Lossless round trip changed the source:\n${JSON.stringify(serialized)}`);
        }

        const types = document.nodes.map(node => node.type);
        if (!types.includes("Comment") || !types.includes("DataRecord")
            || !types.includes("Command") || !types.includes("Unknown")) {
          throw new Error(`Expected typed ordered AST nodes, got ${JSON.stringify(types)}`);
        }
        const engForm = document.nodes.find(node => node.type === "DataRecord" && node.label === "EngForm");
        if (!engForm || engForm.continuationRows.length !== 4) {
          throw new Error(`Expected the EngForm block to own four continuation rows: ${JSON.stringify(engForm)}`);
        }
        if (!engForm.sourceSpan || engForm.sourceSpan.end <= engForm.sourceSpan.start) {
          throw new Error(`Expected a source span on EngForm: ${JSON.stringify(engForm)}`);
        }
        '''
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_lossless_patch_changes_only_last_effective_repeated_record():
    script = textwrap.dedent(
        r'''
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
          "ScatContSym 'SG'  # consumed by the first Scat\r\n",
          "Scat\r\n",
          "UnknownRecord   keep spacing\r\n",
          "ScatContSym 'SU'  # effective table value\r\n",
          "Scat\r\n",
        ].join("");
        const patched = patchEPolyScatInputScript(
          original,
          { scatContSym: "B2U" },
          [],
          { changedKeys: ["scatContSym"] },
        );

        if (!patched.includes("ScatContSym 'SG'  # consumed by the first Scat\r\n")) {
          throw new Error(`The earlier occurrence changed:\n${patched}`);
        }
        if (!patched.includes("ScatContSym 'B2U'")) {
          throw new Error(`The last occurrence was not updated:\n${patched}`);
        }
        if (patched.includes("ScatContSym 'SU'")) {
          throw new Error(`The old effective value remains:\n${patched}`);
        }
        if (!patched.includes("UnknownRecord   keep spacing\r\n")) {
          throw new Error(`Unknown source text changed:\n${patched}`);
        }
        '''
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_ordered_sequence_editor_mutations_preserve_nonsemantic_source_and_order():
    script = textwrap.dedent(
        r'''
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
          appendEPolyScatSequenceNode,
          listEPolyScatSequenceNodes,
          moveEPolyScatSequenceNode,
          removeEPolyScatSequenceNode,
          replaceEPolyScatSequenceNode,
        } = moduleObject.exports;
        const original = [
          "# ordered calculation\r\n",
          "ScatContSym 'SG'  # first symmetry\r\n",
          "Scat\r\n",
          "# keep this separator exactly\r\n",
          "FutureRecord   keep spacing\r\n",
          "ScatContSym 'SU'\r\n",
          "Scat",
        ].join("");

        const sequence = listEPolyScatSequenceNodes(original);
        const labels = sequence.map(item => `${item.type}:${item.label}:${item.occurrence}`);
        const expectedLabels = [
          "DataRecord:ScatContSym:1",
          "Command:Scat:1",
          "Unknown:FutureRecord:1",
          "DataRecord:ScatContSym:2",
          "Command:Scat:2",
        ];
        if (JSON.stringify(labels) !== JSON.stringify(expectedLabels)) {
          throw new Error(`Unexpected ordered projection: ${JSON.stringify(sequence, null, 2)}`);
        }

        const secondSymmetry = sequence[3];
        const replaced = replaceEPolyScatSequenceNode(
          original,
          secondSymmetry.nodeIndex,
          "ScatContSym 'B2U'  # replacement",
        );
        if (!replaced.includes("ScatContSym 'SG'  # first symmetry\r\n")) {
          throw new Error(`Earlier repeated record changed:\n${replaced}`);
        }
        if (!replaced.includes("ScatContSym 'B2U'  # replacement\r\n")) {
          throw new Error(`Replacement did not retain CRLF position:\n${replaced}`);
        }

        const moved = moveEPolyScatSequenceNode(original, sequence[4].nodeIndex, "up");
        const expectedMovedFragment = [
          "# keep this separator exactly\r\n",
          "FutureRecord   keep spacing\r\n",
          "Scat\r\n",
          "ScatContSym 'SU'",
        ].join("");
        if (!moved.includes(expectedMovedFragment)) {
          throw new Error(`Move changed separators or terminal line endings:\n${moved}`);
        }

        const removed = removeEPolyScatSequenceNode(original, sequence[2].nodeIndex);
        if (removed.includes("FutureRecord")) {
          throw new Error(`Unknown record was not removed:\n${removed}`);
        }
        if (!removed.includes("# keep this separator exactly\r\nScatContSym 'SU'")) {
          throw new Error(`Comments around removed record changed:\n${removed}`);
        }

        const appended = appendEPolyScatSequenceNode(original, {
          type: "Command",
          label: "TotalCrossSection",
        });
        if (!appended.endsWith("Scat\r\nTotalCrossSection")) {
          throw new Error(`Append did not follow the document newline convention:\n${appended}`);
        }
        '''
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_structured_sequence_editor_projects_arguments_and_semantic_rows_losslessly():
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
          listEPolyScatSequenceNodes,
          removeEPolyScatSequenceContinuationRow,
          updateEPolyScatSequenceContinuationRow,
          updateEPolyScatSequenceNodeArguments,
        } = moduleObject.exports;
        const original = [
          "Convert '$pt/test03.molden2012' 'molden' # keep convert note",
          "EngForm",
          " 0 2 # charge and formula type",
          " # keep block note",
          " 2",
          " 2.0 -1.0 0",
          " 2.0 -0.5 1",
          "Scat",
          "",
        ].join("\\n");
        const sequence = listEPolyScatSequenceNodes(original);
        const convert = sequence.find(item => item.label === "Convert");
        const engForm = sequence.find(item => item.label === "EngForm");

        if (convert.argumentsText !== "'$pt/test03.molden2012' 'molden'") {
          throw new Error(`Quoted arguments were not preserved: ${convert.argumentsText}`);
        }
        if (JSON.stringify(engForm.continuationRows.map(row => row.value)) !== JSON.stringify([
          "0 2", "2", "2.0 -1.0 0", "2.0 -0.5 1",
        ])) {
          throw new Error(`Unexpected semantic rows: ${JSON.stringify(engForm.continuationRows)}`);
        }
        if (JSON.stringify(engForm.continuationRows.map(row => row.removable)) !== JSON.stringify([
          false, false, true, true,
        ])) {
          throw new Error(`EngForm structural rows were not protected: ${JSON.stringify(engForm.continuationRows)}`);
        }
        try {
          removeEPolyScatSequenceContinuationRow(
            original,
            engForm.nodeIndex,
            engForm.continuationRows[1].sourceLineIndex,
          );
          throw new Error("Removing the EngForm term-count row should fail.");
        } catch (error) {
          if (!String(error.message).includes("structural continuation row")) {
            throw error;
          }
        }

        const updatedConvert = updateEPolyScatSequenceNodeArguments(
          original,
          convert.nodeIndex,
          "'$pt/replacement.molden' 'molden'",
        );
        if (!updatedConvert.includes(
          "Convert '$pt/replacement.molden' 'molden' # keep convert note"
        )) {
          throw new Error(`Header update lost source context:\n${updatedConvert}`);
        }

        const reparsedEngForm = listEPolyScatSequenceNodes(updatedConvert)
          .find(item => item.label === "EngForm");
        const updatedRow = updateEPolyScatSequenceContinuationRow(
          updatedConvert,
          reparsedEngForm.nodeIndex,
          reparsedEngForm.continuationRows[0].sourceLineIndex,
          "1 2",
        );
        if (!updatedRow.includes(" 1 2 # charge and formula type")) {
          throw new Error(`Continuation update lost its inline comment:\n${updatedRow}`);
        }
        if (!updatedRow.includes(" # keep block note")) {
          throw new Error(`Continuation update removed a comment line:\n${updatedRow}`);
        }
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_structured_sequence_editor_adds_and_removes_rows_without_rewriting_neighbors():
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
          appendEPolyScatSequenceContinuationRow,
          listEPolyScatSequenceNodes,
          removeEPolyScatSequenceContinuationRow,
        } = moduleObject.exports;
        const original = [
          "EngForm",
          " 0 1",
          " # preserve between rows",
          " 1",
          " 2.0 -1.0",
          "CustomLegacyRecord keep me",
          "",
          "Scat",
          "",
        ].join("\\n");
        const engForm = listEPolyScatSequenceNodes(original)
          .find(item => item.label === "EngForm");
        const appended = appendEPolyScatSequenceContinuationRow(
          original,
          engForm.nodeIndex,
          "3.0 -0.5",
        );
        const appendedEngForm = listEPolyScatSequenceNodes(appended)
          .find(item => item.label === "EngForm");

        if (!appended.includes(" 2.0 -1.0\\n 3.0 -0.5\\nCustomLegacyRecord keep me")) {
          throw new Error(`New row was not inserted at the end of the block:\n${appended}`);
        }

        const addedRow = appendedEngForm.continuationRows
          .find(row => row.value === "3.0 -0.5");
        const removed = removeEPolyScatSequenceContinuationRow(
          appended,
          appendedEngForm.nodeIndex,
          addedRow.sourceLineIndex,
        );
        if (removed !== original) {
          throw new Error(`Add/remove round trip rewrote neighboring source:\n${removed}`);
        }
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_canonical_schema_covers_all_manual_sample_data_records():
    script = textwrap.dedent(
        r'''
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

        const { EPOLYSCAT_INPUT_SCHEMA } = moduleObject.exports;
        const expected = `AsyPol CnvOrbSel DPotEng EMax EngForm EpsAsym FegeEng FixCharge GrnType
          IPot InitSpinDeg InitSym IterMax JMax JPPMax LFTimeDelayAngles LMax LMaxA
          LMaxI LMaxK Label MFTimeDelayAngles NECenter OrbOcc OrbOccInit PCutRd
          PlaneWvCharge PosFile PosFitL PosGridTol PosPlot PrintBlm ResSearchEng
          RotConstants ScatContSym ScatEng ScatEngN ScatSym SpinDeg TargSpinDeg
          TargSym TestOut VCorr VibAveNInp ViewOrbGrid`.split(/\s+/);
        const labels = new Set(EPOLYSCAT_INPUT_SCHEMA.dataRecords.map(record => record.label));
        const missing = expected.filter(label => !labels.has(label));
        if (missing.length) {
          throw new Error(`Missing manual data records: ${missing.join(", ")}`);
        }
        '''
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_legacy_table_adapter_reads_and_writes_through_lossless_document_model():
    script = textwrap.dedent(
        r'''
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

        const { createLegacyEPolyScatTableObject } = moduleObject.exports;
        let contents = [
          "ScatContSym 'SG' # first state",
          "Scat",
          "ScatContSym 'SU' # table projects the effective state",
          "Scat",
          "FutureRecord untouched",
        ].join("\n");
        const tableObject = createLegacyEPolyScatTableObject({
          getContents: async () => contents,
          setContents: async value => { contents = value; },
        });
        const statePage = tableObject.pages.find(page => page.name === "State Definitions");
        const symmetry = statePage.data.find(field => field.name === "ScatContSym");

        (async () => {
          const before = await symmetry.get();
          if (before !== "SU") {
            throw new Error(`Expected the effective value SU, got ${before}`);
          }
          await symmetry.set("B2U");
          if (!contents.includes("ScatContSym 'SG' # first state")) {
            throw new Error(`Earlier state was changed:\n${contents}`);
          }
          if (!contents.includes("ScatContSym 'B2U'")) {
            throw new Error(`Effective state was not changed:\n${contents}`);
          }
          if (!contents.includes("FutureRecord untouched")) {
            throw new Error(`Unknown source was not preserved:\n${contents}`);
          }
        })().catch(error => {
          console.error(error);
          process.exitCode = 1;
        });
        '''
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_lossless_round_trip_for_all_cached_official_sample_inputs():
    manual_text_dir = (
        APP_DIR.parent.parent
        / "docs"
        / "rrlucchese-epolyscat-manual"
        / "text"
    )
    sample_files = sorted(manual_text_dir.glob("Tests_test??_inp-*.txt"))
    if not sample_files:
        pytest.skip("The local ePolyScat manual cache is not available.")
    assert len(sample_files) == 38

    script = textwrap.dedent(
        f'''
        const fs = require("fs");
        const path = require("path");
        const babel = require("@babel/core");

        const sourcePath = path.join(process.cwd(), "src", "utils", "epolyscat-input-script.js");
        const source = fs.readFileSync(sourcePath, "utf8");
        const compiled = babel.transformSync(source, {{
          plugins: ["@babel/plugin-transform-modules-commonjs"],
        }}).code;
        const moduleObject = {{ exports: {{}} }};
        new Function("module", "exports", compiled)(moduleObject, moduleObject.exports);
        const {{ parseEPolyScatDocument, serializeEPolyScatDocument }} = moduleObject.exports;
        const files = {repr([str(path) for path in sample_files])};

        files.forEach(filename => {{
          const original = fs.readFileSync(filename, "utf8");
          const document = parseEPolyScatDocument(original);
          const serialized = serializeEPolyScatDocument(document);
          if (serialized !== original) {{
            throw new Error(`Round trip changed ${{path.basename(filename)}}`);
          }}
          if (!document.nodes.some(node => node.type === "Command")) {{
            throw new Error(`No commands were recognized in ${{path.basename(filename)}}`);
          }}
          const unknown = document.nodes.filter(node => (
            node.type === "Unknown"
            && !/^test\\d+\\.inp$/i.test(node.sourceLines[0].raw.trim())
          ));
          if (unknown.length) {{
            throw new Error(
              `Unclassified statements in ${{path.basename(filename)}}: `
              + unknown.map(node => node.sourceLines[0].raw).join(" | ")
            );
          }}
        }});
        '''
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_generated_type_zero_engform_and_output_order_follow_manual_semantics():
    script = textwrap.dedent(
        r'''
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
        const { buildEngFormBlock, buildEPolyScatInputScript } = moduleObject.exports;

        const engForm = buildEngFormBlock({
          engFormCharge: 0,
          engFormType: 0,
          engFormTerms: 3,
        });
        if (engForm !== "EngForm\n 0 0") {
          throw new Error(`Type 0 EngForm must not contain a term count:\n${engForm}`);
        }

        const generated = buildEPolyScatInputScript({
          title: "test",
          calculationKind: "scattering",
          lMax: 15,
          eMax: 50,
          fegeEng: 13,
          scatEng: "3.0 4.0",
          engFormCharge: 0,
          engFormType: 0,
          engFormTerms: 0,
          vCorr: "PZ",
          asyPolSwitchD: 0.15,
          asyPolTerms: 1,
          asyPolCenter: 1,
          asyPolValue: 17.5,
          scatContSym: "SG",
          lMaxK: 4,
          convertSource: "$pt/input.molden",
          convertFormat: "molden",
          plotDataFile: "cross-sections.dat",
        }, [{
          fileType: "PlotData",
          valueKey: "plotDataFile",
          disposition: "REWIND",
        }]);
        const fileNameIndex = generated.indexOf("FileName 'PlotData'");
        const producerIndex = generated.indexOf("TotalCrossSection");
        if (fileNameIndex < 0 || producerIndex < 0 || fileNameIndex > producerIndex) {
          throw new Error(`PlotData must be declared before TotalCrossSection:\n${generated}`);
        }
        '''
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_engform_term_count_edit_preserves_retained_rows_and_removes_only_excess_rows():
    script = textwrap.dedent(
        r'''
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
          "EngForm",
          " 0 1 # charge and type",
          " 3 # term count",
          " 2.0 -1.0 # first custom term",
          " 1.5 0.25 # second custom term",
          " 0.5 0.75 # remove this term",
          "FutureRecord remains",
        ].join("\n");
        const patched = patchEPolyScatInputScript(
          original,
          { engFormTerms: 2 },
          [],
          { changedKeys: ["engFormTerms"] },
        );

        if (!patched.includes("2 # term count")) {
          throw new Error(`Term count was not updated:\n${patched}`);
        }
        if (!patched.includes("2.0 -1.0 # first custom term")
            || !patched.includes("1.5 0.25 # second custom term")) {
          throw new Error(`Retained custom terms changed:\n${patched}`);
        }
        if (patched.includes("remove this term")) {
          throw new Error(`Excess term was not removed:\n${patched}`);
        }
        if (!patched.includes("FutureRecord remains")) {
          throw new Error(`Following unknown source changed:\n${patched}`);
        }
        '''
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_asypol_edit_targets_semantic_value_row_without_overwriting_comments():
    script = textwrap.dedent(
        r'''
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
        const { patchEPolyScatInputScript, parseEPolyScatInputScript } = moduleObject.exports;
        const original = [
          "AsyPol",
          " 0.15 # switching distance",
          " 1 # number of terms",
          " 1 # center",
          " 1 # spherical term",
          " # type 2 would read the full tensor",
          " 17.50 # polarizability",
          " 3 # crossing rule",
          " # nearest approach",
          " 0 # matching line",
          " # l = 0 radial function",
          "FegeEng 13.0",
        ].join("\n");
        const patched = patchEPolyScatInputScript(
          original,
          { asyPolValue: 30.4 },
          [],
          { changedKeys: ["asyPolValue"] },
        );
        const parsed = parseEPolyScatInputScript(patched);

        if (parsed.asyPolValue !== 30.4) {
          throw new Error(`The semantic polarizability row was not updated:\n${patched}`);
        }
        if (!patched.includes("# type 2 would read the full tensor")
            || !patched.includes("# nearest approach")
            || !patched.includes("# l = 0 radial function")) {
          throw new Error(`An AsyPol comment was overwritten:\n${patched}`);
        }
        if (!patched.includes("30.4 # polarizability")) {
          throw new Error(`The existing value-row comment was not retained:\n${patched}`);
        }

        const tensorOriginal = [
          "AsyPol",
          " 0.15 # switching distance",
          " 1 # number of terms",
          " 0 # molecular center",
          " 0.0 0.0 0.0 # center coordinates",
          " 2 # tensor term",
          " 8.664 8.664 17.904 0.0 0.0 0.0 # full tensor",
          " 3",
          " 0",
          "FegeEng 13.0",
        ].join("\n");
        const tensorPatched = patchEPolyScatInputScript(
          tensorOriginal,
          {
            asyPolSwitchD: 0.2,
            asyPolTerms: 1,
            asyPolCenter: 0,
            asyPolValue: 8.664,
          },
          [],
          { changedKeys: ["asyPolSwitchD"] },
        );
        if (!tensorPatched.includes("8.664 8.664 17.904 0.0 0.0 0.0 # full tensor")) {
          throw new Error(`Changing SwitchD damaged the tensor row:\n${tensorPatched}`);
        }
        if (parseEPolyScatInputScript(tensorPatched).asyPolValue !== 8.664) {
          throw new Error(`The tensor-form value row was not recognized:\n${tensorPatched}`);
        }
        '''
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_engform_type_transition_removes_and_recreates_term_rows_without_touching_following_records():
    script = textwrap.dedent(
        r'''
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
        const { patchEPolyScatInputScript, parseEPolyScatInputScript } = moduleObject.exports;
        const original = [
          "EngForm # energy formulas",
          " 0 1 # charge and type",
          " 2 # term count",
          " 2.0 -1.0 # custom first term",
          " 1.5 0.25 # custom second term",
          "FegeEng 13.0 # following record",
        ].join("\n");
        const typeZero = patchEPolyScatInputScript(
          original,
          { engFormCharge: 0, engFormType: 0, engFormTerms: 2 },
          [],
          { changedKeys: ["engFormType"] },
        );
        const parsedZero = parseEPolyScatInputScript(typeZero);
        if (parsedZero.engFormType !== 0 || parsedZero.engFormTerms !== undefined) {
          throw new Error(`Type 0 retained formula terms:\n${typeZero}`);
        }
        if (typeZero.includes("custom first term") || typeZero.includes("custom second term")) {
          throw new Error(`Type 0 retained obsolete term rows:\n${typeZero}`);
        }
        if (!typeZero.includes("FegeEng 13.0 # following record")) {
          throw new Error(`The following record changed:\n${typeZero}`);
        }

        const typeTwo = patchEPolyScatInputScript(
          typeZero,
          { engFormCharge: 0, engFormType: 2, engFormTerms: 2 },
          [],
          { changedKeys: ["engFormType"] },
        );
        const parsedTwo = parseEPolyScatInputScript(typeTwo);
        if (parsedTwo.engFormType !== 2 || parsedTwo.engFormTerms !== 2) {
          throw new Error(`Type 2 did not recreate its term block:\n${typeTwo}`);
        }
        const defaultTerms = typeTwo.match(/2\.0 -1\.0 0/g) || [];
        if (defaultTerms.length !== 2) {
          throw new Error(`Type 2 expected two default term rows:\n${typeTwo}`);
        }
        if (!typeTwo.includes("FegeEng 13.0 # following record")) {
          throw new Error(`The following record changed after recreation:\n${typeTwo}`);
        }
        '''
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


def test_generated_output_destinations_precede_the_commands_that_write_them():
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

        const { buildEPolyScatInputScript } = moduleObject.exports;
        const baseValues = {
          title: "manual-backed output ordering",
          lMax: 15,
          lMaxI: 15,
          eMax: 50,
          fegeEng: 13,
          scatEng: "3.0 4.0",
          engFormCharge: 0,
          engFormType: 0,
          engFormTerms: 0,
          vCorr: "PZ",
          asyPolSwitchD: 1,
          asyPolTerms: 1,
          asyPolCenter: 1,
          asyPolValue: 10,
          scatContSym: "SG",
          scatSym: "A1",
          lMaxK: 4,
          convertSource: "$pt/input.molden",
          convertFormat: "molden",
          initSym: "A1",
          initSpinDeg: 1,
          orbOccInit: "2 2",
          orbOcc: "2 1",
          spinDeg: 1,
          targSym: "A1",
          targSpinDeg: 1,
          iPot: 1,
          plotDataFile: "cross-sections.dat",
        };
        const outputDefinitions = [{
          fileType: "PlotData",
          valueKey: "plotDataFile",
          disposition: "REWIND",
        }];

        function assertBefore(contents, first, second) {
          const firstIndex = contents.indexOf(first);
          const secondIndex = contents.indexOf(second);
          if (firstIndex < 0 || secondIndex < 0 || firstIndex >= secondIndex) {
            throw new Error(`Expected ${first} before ${second}\n${contents}`);
          }
        }

        const scattering = buildEPolyScatInputScript({
          ...baseValues,
          calculationKind: "scattering",
        }, outputDefinitions);
        assertBefore(scattering, "FileName 'PlotData'", "TotalCrossSection");

        const photoionization = buildEPolyScatInputScript({
          ...baseValues,
          calculationKind: "photoionization",
        }, outputDefinitions);
        assertBefore(photoionization, "FileName 'PlotData'", "GetCro");
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr
