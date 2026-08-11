function normalizedName(file) {
  return String(file && file.name ? file.name : "").toLowerCase();
}

function normalizedFileType(file) {
  return String(file && (file.fileType || file.file_type || file.type) ? (file.fileType || file.file_type || file.type) : "").toLowerCase();
}

function stagedInputReference(file) {
  const fileName = String(file && file.name ? file.name : "");
  const baseName = fileName.split(/[\\/]/).pop();
  return baseName ? `$pt/${baseName}` : "";
}

function findPreferredFile(outputFiles, predicates) {
  const files = Array.isArray(outputFiles) ? outputFiles : [];

  for (const predicate of predicates) {
    const file = files.find(predicate);
    if (file) {
      return file;
    }
  }

  return null;
}

function isGaussianConvertOutput(file) {
  const name = normalizedName(file);

  return name === "gaussian.log"
      || name === "fort.7"
      || /\.g0[39]$/.test(name)
      || /\.g16$/.test(name)
      || (name.includes("gaussian") && name.endsWith(".log"));
}

function isMoldenConvertOutput(file) {
  const name = normalizedName(file);

  return name === "molden.dat"
      || name.includes("molden")
      || name.endsWith(".molden")
      || name.endsWith(".molden.dat");
}

function isDumpOutFile(file) {
  const name = normalizedName(file);
  const fileType = normalizedFileType(file);

  return fileType === "dumpout"
      || name.includes("dumpidy")
      || name.includes("dumpout");
}

function isBendOrientOutput(file) {
  const name = normalizedName(file);
  const fileType = normalizedFileType(file);

  return fileType === "bendorient"
      || name.includes("bendorient");
}

function isOrientNCroOutput(file) {
  const name = normalizedName(file);
  const fileType = normalizedFileType(file);

  return fileType === "orientncro"
      || name.includes("orientncro")
      || name.includes("orinetncro");
}

function isCubeOutput(file) {
  const name = normalizedName(file);

  return name.endsWith(".cube") || name.endsWith(".cub");
}

export function buildEPolyScatInputDataSelectionValues(file, sourceApplicationId = "") {
  const sourceApplication = String(sourceApplicationId || "");
  const useGaussianFormat = sourceApplication === "Gaussian16"
      || (!sourceApplication && isGaussianConvertOutput(file));

  return {
    convertSource: stagedInputReference(file),
    convertFormat: useGaussianFormat ? "gaussian" : "molden",
  };
}

function buildInputBinding(inputFileName, outputFile, dataEntryValues = null) {
  if (!inputFileName || !outputFile) {
    return null;
  }

  return {
    inputFileName,
    outputFile,
    dataEntryValues,
  };
}

function buildEPolyScatInputDataBinding(outputFiles, sourceApplicationId) {
  const sourceApplication = String(sourceApplicationId || "");
  const sourceIsGaussian = sourceApplication === "Gaussian16";
  const sourceIsOpenMolcas = sourceApplication === "OpenMolcas";
  const outputFile = sourceIsGaussian
      ? findPreferredFile(outputFiles, [isGaussianConvertOutput])
      : sourceIsOpenMolcas
          ? findPreferredFile(outputFiles, [isMoldenConvertOutput])
          : findPreferredFile(outputFiles, [isMoldenConvertOutput, isGaussianConvertOutput]);

  if (!outputFile) {
    return null;
  }

  return buildInputBinding(
      "ePolyScat_Input_Data",
      outputFile,
      buildEPolyScatInputDataSelectionValues(outputFile, sourceApplication),
  );
}

function buildAnalysisBinding(outputFiles, targetApplicationId, requiredFileName) {
  const targetApplication = String(targetApplicationId || "");
  const inputFileName = requiredFileName || "";
  let outputFile = null;

  if (targetApplication === "CnvLinFull") {
    outputFile = findPreferredFile(outputFiles, [isDumpOutFile]);
  } else if (targetApplication === "CnvMath" || targetApplication === "CnvMatLab") {
    outputFile = findPreferredFile(outputFiles, [isBendOrientOutput]);
  } else if (targetApplication === "NRFPAD") {
    outputFile = findPreferredFile(outputFiles, [isOrientNCroOutput]);
  } else if (targetApplication === "Cube2igor") {
    outputFile = findPreferredFile(outputFiles, [isCubeOutput]);
  } else if (targetApplication === "MoldenMerge") {
    outputFile = findPreferredFile(outputFiles, [isMoldenConvertOutput]);
  }

  return buildInputBinding(inputFileName, outputFile);
}

function findExactRequiredFile(outputFiles, requiredFileName) {
  const requiredName = String(requiredFileName || "").toLowerCase();
  if (!requiredName) {
    return null;
  }

  return findPreferredFile(outputFiles, [
    file => normalizedName(file) === requiredName,
  ]);
}

export function buildWorkflowOutputInputBinding({
  outputFiles = [],
  targetStageId = "",
  sourceApplicationId = "",
  targetApplicationId = "",
  requiredFileName = "",
} = {}) {
  if (targetStageId === "ePolyScat_Run") {
    return buildEPolyScatInputDataBinding(outputFiles, sourceApplicationId);
  }

  if (targetStageId === "Analysis") {
    const binding = buildAnalysisBinding(outputFiles, targetApplicationId, requiredFileName);
    if (binding) {
      return binding;
    }
  }

  const exactFile = findExactRequiredFile(outputFiles, requiredFileName);
  return buildInputBinding(requiredFileName, exactFile);
}
