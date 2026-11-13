import { patchEPolyScatInputScript } from "./epolyscat-input-script";
import { buildEPolyScatInputDataSelectionValues } from "./workflow-file-linking";

const SCRIPT_INPUT_NAME = "ePolyscat_Input_File";
const RUNTIME_INPUT_NAME = "ePolyScat_Input_Data";

function tutorialFileEntries(sourceRun) {
  return (sourceRun.inputs || [])
      .filter(input => input.type === "files")
      .reduce((entries, input) => entries.concat(
          (input.files || []).map(file => ({ inputFileName: input.name, file })),
      ), []);
}

function localHydratedFile(file, contents) {
  return {
    ...file,
    id: null,
    contents,
    dataProductURI: null,
    dataProductUri: null,
    data_product_uri: null,
    "data-product-uri": null,
    downloadURL: null,
    download_url: null,
    deleted: false,
    generatedByNewRun: false,
    unchanged: false,
    replaceCurrent: true,
  };
}

export function requiredFilesAreReady(requiredFiles = [], inputFiles = {}) {
  return requiredFiles.every(requiredFile => {
    const inputFile = inputFiles[requiredFile.name];
    const uploadedFiles = inputFile && Array.isArray(inputFile.files)
        ? inputFile.files.filter(file => !file.deleted && !file.generatedByNewRun)
        : [];
    return uploadedFiles.length >= (requiredFile.minimumFileCount || 1);
  });
}

export async function prepareTutorialCloneInputs({
  sourceRun,
  fetchContent,
  dataEntryValues = {},
  outputDefinitions = [],
} = {}) {
  const metadata = sourceRun && sourceRun.workflowMetadata
      ? sourceRun.workflowMetadata
      : {};
  const runtimeDataName = metadata.tutorialRuntimeDataName;
  const convertFormat = metadata.tutorialConvertFormat;
  if (!runtimeDataName || !convertFormat) {
    throw new Error("Tutorial input metadata is incomplete.");
  }
  if (typeof fetchContent !== "function") {
    throw new Error("Tutorial input loader is unavailable.");
  }

  const fileEntries = tutorialFileEntries(sourceRun);
  const scriptEntries = fileEntries.filter(
      entry => entry.inputFileName === SCRIPT_INPUT_NAME,
  );
  const runtimeEntries = fileEntries.filter(
      entry => entry.inputFileName === RUNTIME_INPUT_NAME
          && entry.file.name === runtimeDataName,
  );
  if (scriptEntries.length !== 1 || runtimeEntries.length !== 1) {
    throw new Error("Tutorial runtime inputs are incomplete.");
  }

  const hydratedEntries = await Promise.all(
      [scriptEntries[0], runtimeEntries[0]].map(async entry => ({
        ...entry,
        contents: await fetchContent({
          runId: sourceRun.id,
          filename: entry.file.name,
        }),
      })),
  );
  const hydratedScriptEntry = hydratedEntries[0];
  const convertRecords = String(hydratedScriptEntry.contents || "")
      .match(/^\s*Convert\b.*$/gmi) || [];
  if (convertRecords.length !== 1) {
    throw new Error("Tutorial script must contain exactly one Convert command.");
  }

  const linkedValues = buildEPolyScatInputDataSelectionValues(
      runtimeEntries[0].file,
      "",
      convertFormat,
  );
  const scriptContent = patchEPolyScatInputScript(
      hydratedScriptEntry.contents,
      { ...dataEntryValues, ...linkedValues },
      outputDefinitions,
      { changedKeys: ["convertSource", "convertFormat"] },
  );

  return {
    dataEntryValues: linkedValues,
    scriptContent,
    files: hydratedEntries.map(entry => ({
      inputFileName: entry.inputFileName,
      file: localHydratedFile(
          entry.file,
          entry === hydratedScriptEntry ? scriptContent : entry.contents,
      ),
    })),
  };
}
