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


def test_tutorial_clone_hydration_is_atomic_retriable_and_clears_remote_uris():
    script = textwrap.dedent(
        r'''
        const fs = require("fs");
        const path = require("path");
        const babel = require("@babel/core");

        function compileModule(relativePath, dependencies = {}) {
          const source = fs.readFileSync(path.join(process.cwd(), relativePath), "utf8");
          const compiled = babel.transformSync(source, {
            plugins: ["@babel/plugin-transform-modules-commonjs"],
          }).code;
          const moduleObject = { exports: {} };
          const localRequire = name => {
            return name in dependencies ? dependencies[name] : require(name);
          };
          new Function("module", "exports", "require", compiled)(
            moduleObject,
            moduleObject.exports,
            localRequire,
          );
          return moduleObject.exports;
        }

        const inputScript = compileModule(
          path.join("src", "utils", "epolyscat-input-script.js"),
        );
        const workflowFileLinking = compileModule(
          path.join("src", "utils", "workflow-file-linking.js"),
        );
        const tutorialClone = compileModule(
          path.join("src", "utils", "tutorial-clone.js"),
          {
            "./epolyscat-input-script": inputScript,
            "./workflow-file-linking": workflowFileLinking,
          },
        );
        const {
          prepareTutorialCloneInputs,
          requiredFilesAreReady,
        } = tutorialClone;
        const catalogRoot = path.join(
          process.cwd(),
          "data",
          "epolyscat",
          "manual_tutorials",
        );
        const manifest = JSON.parse(
          fs.readFileSync(path.join(catalogRoot, "manifest.json"), "utf8"),
        );
        const requiredFiles = [
          { name: "ePolyscat_Input_File", minimumFileCount: 1 },
          { name: "ePolyScat_Input_Data", minimumFileCount: 1 },
        ];

        function sourceRunFor(tutorialCase) {
          const remoteFields = {
            dataProductURI: "airavata-dp://must-clear",
            data_product_uri: "airavata-dp://must-clear",
            "data-product-uri": "airavata-dp://must-clear",
            downloadURL: "/remote-download/must-clear",
            unchanged: true,
          };
          return {
            id: 9001,
            workflowMetadata: {
              tutorialRuntimeDataName: tutorialCase.runtime_data,
              tutorialConvertFormat: tutorialCase.convert_format,
            },
            inputs: [
              {
                type: "files",
                name: "ePolyscat_Input_File",
                files: [{ name: tutorialCase.script, ...remoteFields }],
              },
              {
                type: "files",
                name: "ePolyScat_Input_Data",
                files: [{ name: tutorialCase.runtime_data, ...remoteFields }],
              },
            ],
          };
        }

        async function main() {
          for (const tutorialCase of manifest.cases) {
            const fetchCalls = [];
            const result = await prepareTutorialCloneInputs({
              sourceRun: sourceRunFor(tutorialCase),
              fetchContent: async ({ filename }) => {
                fetchCalls.push(filename);
                return fs.readFileSync(
                  path.join(catalogRoot, tutorialCase.id, filename),
                  "utf8",
                );
              },
              dataEntryValues: {},
              outputDefinitions: [],
            });
            const expectedFetches = [tutorialCase.script, tutorialCase.runtime_data].sort();
            if (fetchCalls.slice().sort().join("|") !== expectedFetches.join("|")) {
              throw new Error(`${tutorialCase.id}: hydration did not fetch exactly both inputs`);
            }
            const parsed = inputScript.parseEPolyScatInputScript(result.scriptContent);
            if (parsed.convertSource !== `$pt/${tutorialCase.runtime_data}`) {
              throw new Error(`${tutorialCase.id}: wrong Convert source ${parsed.convertSource}`);
            }
            if (parsed.convertFormat !== tutorialCase.convert_format) {
              throw new Error(`${tutorialCase.id}: wrong Convert format ${parsed.convertFormat}`);
            }
            if (result.files.length !== 2) {
              throw new Error(`${tutorialCase.id}: expected exactly two hydrated inputs`);
            }
            for (const entry of result.files) {
              const file = entry.file;
              for (const field of [
                "dataProductURI",
                "data_product_uri",
                "data-product-uri",
                "downloadURL",
              ]) {
                if (file[field] !== null) {
                  throw new Error(`${tutorialCase.id}: ${field} was not cleared`);
                }
              }
              if (file.unchanged !== false || file.replaceCurrent !== true) {
                throw new Error(`${tutorialCase.id}: hydrated file is not a local replacement`);
              }
            }

            const inputFiles = Object.fromEntries(result.files.map(entry => [
              entry.inputFileName,
              { files: [entry.file] },
            ]));
            if (!requiredFilesAreReady(requiredFiles, inputFiles)) {
              throw new Error(`${tutorialCase.id}: both hydrated files should pass save gating`);
            }
            inputFiles.ePolyScat_Input_Data.files[0].deleted = true;
            if (requiredFilesAreReady(requiredFiles, inputFiles)) {
              throw new Error(`${tutorialCase.id}: removing runtime data must block save`);
            }
          }

          const retryCase = manifest.cases[1];
          let shouldFail = true;
          let committedFiles = [];
          const fetchContent = async ({ filename }) => {
            if (shouldFail && filename === retryCase.runtime_data) {
              throw new Error("simulated runtime fetch failure");
            }
            return fs.readFileSync(
              path.join(catalogRoot, retryCase.id, filename),
              "utf8",
            );
          };
          try {
            const firstResult = await prepareTutorialCloneInputs({
              sourceRun: sourceRunFor(retryCase),
              fetchContent,
            });
            committedFiles = firstResult.files;
            throw new Error("expected initial hydration to fail");
          } catch (error) {
            if (!String(error.message).includes("simulated runtime fetch failure")) {
              throw error;
            }
          }
          if (committedFiles.length !== 0) {
            throw new Error("failed atomic hydration exposed a partial file result");
          }
          shouldFail = false;
          const retryResult = await prepareTutorialCloneInputs({
            sourceRun: sourceRunFor(retryCase),
            fetchContent,
          });
          if (retryResult.files.length !== 2) {
            throw new Error("retry did not hydrate both inputs");
          }
        }

        main().catch(error => {
          console.error(error);
          process.exitCode = 1;
        });
        '''
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr
