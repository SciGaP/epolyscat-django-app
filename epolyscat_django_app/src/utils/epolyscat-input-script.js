const FIELD_RECORDS = {
  title: ["Title"],
  lMax: ["LMax"],
  lMaxI: ["LMaxI"],
  eMax: ["EMax"],
  fegeEng: ["FegeEng"],
  scatEng: ["ScatEng", "Scat"],
  lMaxK: ["LMaxK"],
  vCorr: ["VCorr"],
  scatContSym: ["ScatContSym"],
  scatSym: ["ScatSym"],
  initSym: ["InitSym"],
  initSpinDeg: ["InitSpinDeg"],
  orbOccInit: ["OrbOccInit"],
  orbOcc: ["OrbOcc"],
  spinDeg: ["SpinDeg"],
  targSym: ["TargSym"],
  targSpinDeg: ["TargSpinDeg"],
  iPot: ["IPot"],
  convertSource: ["Convert"],
  convertFormat: ["Convert"],
  engFormCharge: ["EngForm"],
  engFormType: ["EngForm"],
  engFormTerms: ["EngForm"],
  asyPolSwitchD: ["AsyPol"],
  asyPolTerms: ["AsyPol"],
  asyPolCenter: ["AsyPol"],
  asyPolValue: ["AsyPol"],
};

const SINGLE_LINE_RECORD_ORDER = [
  "LMax",
  "LMaxI",
  "EMax",
  "FegeEng",
  "ScatEng",
  "InitSym",
  "InitSpinDeg",
  "OrbOccInit",
  "OrbOcc",
  "SpinDeg",
  "TargSym",
  "TargSpinDeg",
  "IPot",
  "VCorr",
  "ScatContSym",
  "ScatSym",
  "LMaxK",
  "Convert",
];

function hasValue(value) {
  return value !== undefined && value !== null && value !== "";
}

export function quoteEPolyScatValue(value) {
  return `'${String(value == null ? "" : value).replace(/'/g, "''")}'`;
}

export function stripEPolyScatComment(line) {
  let inQuote = false;

  for (let index = 0; index < line.length; index++) {
    const char = line[index];
    const nextChar = line[index + 1];

    if (char === "'" && inQuote && nextChar === "'") {
      index++;
      continue;
    }

    if (char === "'") {
      inQuote = !inQuote;
      continue;
    }

    if (char === "#" && !inQuote) {
      return line.slice(0, index);
    }
  }

  return line;
}

export function tokenizeEPolyScatLine(line) {
  const withoutComment = stripEPolyScatComment(line);
  const tokens = withoutComment.match(/'(?:[^']|'')*'|\S+/g) || [];

  return tokens.map(token => {
    if (token.startsWith("'") && token.endsWith("'")) {
      return token.slice(1, -1).replace(/''/g, "'");
    }

    return token;
  });
}

export function normalizeEPolyScatRecord(record) {
  const recordMap = {
    lmax: "LMax",
    lmaxi: "LMaxI",
    emax: "EMax",
    fegeeng: "FegeEng",
    scateng: "ScatEng",
    scat: "Scat",
    lmaxk: "LMaxK",
    vcorr: "VCorr",
    scatcontsym: "ScatContSym",
    scatsym: "ScatSym",
    initsym: "InitSym",
    initspindeg: "InitSpinDeg",
    orboccinit: "OrbOccInit",
    orbocc: "OrbOcc",
    spindeg: "SpinDeg",
    targsym: "TargSym",
    targspindeg: "TargSpinDeg",
    ipot: "IPot",
    convert: "Convert",
    filename: "FileName",
    engform: "EngForm",
    asypol: "AsyPol",
    genformphion: "GenFormPhIon",
    dipoleop: "DipoleOp",
    phion: "PhIon",
    getcro: "GetCro",
  };

  return recordMap[String(record || "").toLowerCase()] || record;
}

export function normalizeEPolyScatInputContents(contents) {
  const text = String(contents == null ? "" : contents);
  const normalizedText = /^\s*\{\\rtf/i.test(text) ? extractTextFromRtf(text) : text;
  return normalizeEscapedTextArtifacts(normalizedText);
}

function normalizeEscapedTextArtifacts(contents) {
  return String(contents)
      .replace(/\\(?=\s)/g, "")
      .replace(/\\(?=\r?\n|$)/g, "");
}

function extractTextFromRtf(contents) {
  const contentStart = contents.search(/\\pard\b/);
  let text = contentStart >= 0 ? contents.slice(contentStart) : contents;

  text = text
      .replace(/\\par(?![a-zA-Z])/g, "\n")
      .replace(/\\line(?![a-zA-Z])/g, "\n")
      .replace(/\\tab(?![a-zA-Z])/g, "\t")
      .replace(/\\'([0-9a-fA-F]{2})/g, (match, hex) => String.fromCharCode(parseInt(hex, 16)))
      .replace(/\\u(-?\d+)\??/g, (match, code) => {
        let value = Number(code);
        if (value < 0) {
          value += 65536;
        }
        return String.fromCharCode(value);
      })
      .replace(/\\([{}\\])/g, "$1")
      .replace(/\\ /g, " ")
      .replace(/\\[a-zA-Z]+-?\d* ?/g, "")
      .replace(/\\[^a-zA-Z0-9\s]/g, "")
      .replace(/[{}]/g, "");

  return text
      .split(/\r?\n/)
      .map(line => line.trimEnd())
      .join("\n")
      .replace(/\n{3,}/g, "\n\n")
      .trim();
}

export function parseNumericToken(value) {
  const numericValue = Number(value);
  return value !== "" && value != null && !Number.isNaN(numericValue) ? numericValue : value;
}

function setParsedValue(values, key, value) {
  if (value !== undefined) {
    values[key] = value;
  }
}

function isNumericToken(value) {
  return value !== "" && value != null && !Number.isNaN(Number(value));
}

function collectFollowingValueTokens(lines, startIndex, maxRows) {
  const valueRows = [];
  let nextIndex = startIndex;

  for (let index = startIndex + 1; index < lines.length && valueRows.length < maxRows; index++) {
    const tokens = tokenizeEPolyScatLine(lines[index]);

    if (tokens.length === 0) {
      if (valueRows.length > 0) {
        break;
      }
      continue;
    }

    if (!isNumericToken(tokens[0])) {
      break;
    }

    valueRows.push(tokens);
    nextIndex = index;
  }

  return { valueRows, nextIndex };
}

function parseConvertRecord(tokens) {
  return {
    convertSource: tokens[1] || "",
    convertFormat: tokens[2] || "",
  };
}

function parseFileNameRecord(tokens, outputDefinitions) {
  const fileType = tokens[1];
  const output = outputDefinitions.find(definition => definition.fileType === fileType);

  if (!output) {
    return {};
  }

  return {
    [output.valueKey]: tokens[2] || "",
  };
}

function parseEngFormBlock(lines, startIndex) {
  const records = collectFollowingValueTokens(lines, startIndex, 2);
  const [formulaRecord, termsRecord] = records.valueRows;
  const values = {};

  if (formulaRecord) {
    setParsedValue(values, "engFormCharge", parseNumericToken(formulaRecord[0]));
    setParsedValue(values, "engFormType", parseNumericToken(formulaRecord[1]));
  }

  if (termsRecord) {
    setParsedValue(values, "engFormTerms", parseNumericToken(termsRecord[0]));
  }

  return {
    values,
    nextIndex: records.nextIndex,
  };
}

function parseAsyPolBlock(lines, startIndex) {
  const records = collectFollowingValueTokens(lines, startIndex, 7);
  const [switchRecord, termsRecord, centerRecord, , polarizabilityRecord] = records.valueRows;
  const values = {};

  if (switchRecord) {
    setParsedValue(values, "asyPolSwitchD", parseNumericToken(switchRecord[0]));
  }

  if (termsRecord) {
    setParsedValue(values, "asyPolTerms", parseNumericToken(termsRecord[0]));
  }

  if (centerRecord) {
    setParsedValue(values, "asyPolCenter", parseNumericToken(centerRecord[0]));
  }

  if (polarizabilityRecord) {
    setParsedValue(values, "asyPolValue", parseNumericToken(polarizabilityRecord[0]));
  }

  return {
    values,
    nextIndex: records.nextIndex,
  };
}

export function parseEPolyScatInputScript(contents, outputDefinitions = []) {
  const lines = normalizeEPolyScatInputContents(contents).split(/\r?\n/);
  const parsedValues = {};

  for (let index = 0; index < lines.length; index++) {
    const rawLine = lines[index];
    const trimmedLine = rawLine.trim();

    if (!parsedValues.title && trimmedLine.startsWith("#")) {
      const title = trimmedLine.replace(/^#+\s*/, "").trim();
      if (title) {
        parsedValues.title = title;
      }
    }

    const tokens = tokenizeEPolyScatLine(rawLine);
    if (tokens.length === 0) {
      continue;
    }

    const record = normalizeEPolyScatRecord(tokens[0]);
    if (record === "LMax") {
      parsedValues.lMax = parseNumericToken(tokens[1]);
    } else if (record === "LMaxI") {
      parsedValues.lMaxI = parseNumericToken(tokens[1]);
    } else if (record === "EMax") {
      parsedValues.eMax = parseNumericToken(tokens[1]);
    } else if (record === "FegeEng") {
      parsedValues.fegeEng = parseNumericToken(tokens[1]);
    } else if (record === "ScatEng") {
      parsedValues.scatEng = tokens.slice(1).join(" ");
      parsedValues.calculationKind = "scattering";
    } else if (record === "Scat") {
      parsedValues.calculationKind = "scattering";
      if (tokens.length > 1) {
        parsedValues.scatEng = tokens.slice(1).join(" ");
      }
    } else if (record === "LMaxK") {
      parsedValues.lMaxK = parseNumericToken(tokens[1]);
    } else if (record === "VCorr") {
      parsedValues.vCorr = tokens[1] || "";
    } else if (record === "ScatContSym") {
      parsedValues.scatContSym = tokens[1] || "";
    } else if (record === "ScatSym") {
      parsedValues.scatSym = tokens[1] || "";
    } else if (record === "InitSym") {
      parsedValues.initSym = tokens[1] || "";
      parsedValues.calculationKind = "photoionization";
    } else if (record === "InitSpinDeg") {
      parsedValues.initSpinDeg = parseNumericToken(tokens[1]);
      parsedValues.calculationKind = "photoionization";
    } else if (record === "OrbOccInit") {
      parsedValues.orbOccInit = tokens.slice(1).join(" ");
      parsedValues.calculationKind = "photoionization";
    } else if (record === "OrbOcc") {
      parsedValues.orbOcc = tokens.slice(1).join(" ");
      parsedValues.calculationKind = "photoionization";
    } else if (record === "SpinDeg") {
      parsedValues.spinDeg = parseNumericToken(tokens[1]);
      parsedValues.calculationKind = "photoionization";
    } else if (record === "TargSym") {
      parsedValues.targSym = tokens[1] || "";
      parsedValues.calculationKind = "photoionization";
    } else if (record === "TargSpinDeg") {
      parsedValues.targSpinDeg = parseNumericToken(tokens[1]);
      parsedValues.calculationKind = "photoionization";
    } else if (record === "IPot") {
      parsedValues.iPot = parseNumericToken(tokens[1]);
      parsedValues.calculationKind = "photoionization";
    } else if (record === "Convert") {
      Object.assign(parsedValues, parseConvertRecord(tokens));
    } else if (record === "FileName") {
      Object.assign(parsedValues, parseFileNameRecord(tokens, outputDefinitions));
    } else if (record === "EngForm") {
      const parsedBlock = parseEngFormBlock(lines, index);
      Object.assign(parsedValues, parsedBlock.values);
      index = parsedBlock.nextIndex;
    } else if (record === "AsyPol") {
      const parsedBlock = parseAsyPolBlock(lines, index);
      Object.assign(parsedValues, parsedBlock.values);
      index = parsedBlock.nextIndex;
    } else if (["GenFormPhIon", "DipoleOp", "PhIon", "GetCro"].includes(record)) {
      parsedValues.calculationKind = "photoionization";
    }
  }

  return parsedValues;
}

function buildFileNameRecord(output, values) {
  const fileName = values[output.valueKey];
  if (!hasValue(fileName)) {
    return null;
  }

  return [
    "FileName",
    quoteEPolyScatValue(output.fileType),
    quoteEPolyScatValue(fileName),
    output.disposition ? quoteEPolyScatValue(output.disposition) : "",
  ].filter(Boolean).join(" ");
}

function buildOutputFileNameRecords(values, outputDefinitions, fileTypes = []) {
  return outputDefinitions
      .filter(output => fileTypes.length === 0 || fileTypes.includes(output.fileType))
      .map(output => buildFileNameRecord(output, values))
      .filter(Boolean);
}

export function buildEngFormBlock(values) {
  const termCount = Number(values.engFormTerms) || 0;
  const terms = Array.from({ length: termCount }, () => " 2.0 -1.0");
  return [
    "EngForm",
    ` ${values.engFormCharge} ${values.engFormType}`,
    ` ${termCount}`,
    ...terms,
  ].join("\n");
}

export function buildAsyPolBlock(values) {
  return [
    "AsyPol",
    ` ${values.asyPolSwitchD}`,
    ` ${values.asyPolTerms}`,
    ` ${values.asyPolCenter}`,
    " 1",
    ` ${values.asyPolValue}`,
    " 3",
    " 0",
  ].join("\n");
}

export function buildEPolyScatInputScript(values, outputDefinitions = []) {
  const lines = [
    `# ${values.title}`,
    `LMax ${values.lMax}`,
  ];

  if (hasValue(values.lMaxI)) {
    lines.push(`LMaxI ${values.lMaxI}`);
  }

  lines.push(
      `EMax ${values.eMax}`,
      `FegeEng ${values.fegeEng}`,
      `ScatEng ${values.scatEng}`
  );

  if (values.calculationKind === "photoionization") {
    lines.push(
        `InitSym ${quoteEPolyScatValue(values.initSym)}`,
        `InitSpinDeg ${values.initSpinDeg}`,
        `OrbOccInit ${values.orbOccInit}`,
        `OrbOcc ${values.orbOcc}`,
        `SpinDeg ${values.spinDeg}`,
        `TargSym ${quoteEPolyScatValue(values.targSym)}`,
        `TargSpinDeg ${values.targSpinDeg}`,
        `IPot ${values.iPot}`
    );
  } else {
    lines.push(
        buildEngFormBlock(values),
        `VCorr ${quoteEPolyScatValue(values.vCorr)}`,
        buildAsyPolBlock(values),
        `ScatContSym ${quoteEPolyScatValue(values.scatContSym)}`,
        `LMaxK ${values.lMaxK}`
    );
  }

  lines.push(
      "",
      `Convert ${quoteEPolyScatValue(values.convertSource)} ${quoteEPolyScatValue(values.convertFormat)}`,
      "GetBlms",
      "ExpOrb"
  );

  if (values.calculationKind === "photoionization") {
    lines.push(
        "",
        `ScatSym ${quoteEPolyScatValue(values.scatSym)}`,
        `ScatContSym ${quoteEPolyScatValue(values.scatContSym)}`,
        ...buildOutputFileNameRecords(values, outputDefinitions, ["MatrixElements"]),
        "GenFormPhIon",
        "DipoleOp",
        "GetPot",
        "PhIon",
        ...buildOutputFileNameRecords(values, outputDefinitions, ["PlotData"]),
        "GetCro"
    );
  } else {
    lines.push(
        "GetPot",
        `Scat ${values.scatEng}`,
        ...buildOutputFileNameRecords(values, outputDefinitions, ["PlotData"]),
        "TotalCrossSection"
    );
  }

  return lines.join("\n").replace(/\n{3,}/g, "\n\n");
}

function buildSingleLineReplacements(values) {
  const replacements = {};

  if (hasValue(values.lMax)) replacements.LMax = `LMax ${values.lMax}`;
  if (hasValue(values.lMaxI)) replacements.LMaxI = `LMaxI ${values.lMaxI}`;
  if (hasValue(values.eMax)) replacements.EMax = `EMax ${values.eMax}`;
  if (hasValue(values.fegeEng)) replacements.FegeEng = `FegeEng ${values.fegeEng}`;
  if (hasValue(values.scatEng)) {
    replacements.ScatEng = `ScatEng ${values.scatEng}`;
    replacements.Scat = `Scat ${values.scatEng}`;
  }
  if (hasValue(values.initSym)) replacements.InitSym = `InitSym ${quoteEPolyScatValue(values.initSym)}`;
  if (hasValue(values.initSpinDeg)) replacements.InitSpinDeg = `InitSpinDeg ${values.initSpinDeg}`;
  if (hasValue(values.orbOccInit)) replacements.OrbOccInit = `OrbOccInit ${values.orbOccInit}`;
  if (hasValue(values.orbOcc)) replacements.OrbOcc = `OrbOcc ${values.orbOcc}`;
  if (hasValue(values.spinDeg)) replacements.SpinDeg = `SpinDeg ${values.spinDeg}`;
  if (hasValue(values.targSym)) replacements.TargSym = `TargSym ${quoteEPolyScatValue(values.targSym)}`;
  if (hasValue(values.targSpinDeg)) replacements.TargSpinDeg = `TargSpinDeg ${values.targSpinDeg}`;
  if (hasValue(values.iPot)) replacements.IPot = `IPot ${values.iPot}`;
  if (hasValue(values.vCorr)) replacements.VCorr = `VCorr ${quoteEPolyScatValue(values.vCorr)}`;
  if (hasValue(values.scatContSym)) replacements.ScatContSym = `ScatContSym ${quoteEPolyScatValue(values.scatContSym)}`;
  if (hasValue(values.scatSym)) replacements.ScatSym = `ScatSym ${quoteEPolyScatValue(values.scatSym)}`;
  if (hasValue(values.lMaxK)) replacements.LMaxK = `LMaxK ${values.lMaxK}`;
  if (hasValue(values.convertSource) || hasValue(values.convertFormat)) {
    replacements.Convert = `Convert ${quoteEPolyScatValue(values.convertSource)} ${quoteEPolyScatValue(values.convertFormat)}`;
  }

  return replacements;
}

function buildFileNameReplacements(values, outputDefinitions) {
  return outputDefinitions.reduce((replacements, output) => {
    const record = buildFileNameRecord(output, values);
    if (record) {
      replacements[output.fileType] = record;
    }

    return replacements;
  }, {});
}

function getTargetPatch(options, outputDefinitions) {
  const changedKeys = options && Array.isArray(options.changedKeys) ? options.changedKeys : [];
  if (changedKeys.length === 0) {
    return null;
  }

  const records = new Set();
  const fileTypes = new Set();
  changedKeys.forEach(key => {
    const mappedRecords = FIELD_RECORDS[key] || [];
    mappedRecords.forEach(record => records.add(record));
    outputDefinitions
        .filter(output => output.valueKey === key)
        .forEach(output => fileTypes.add(output.fileType));
  });

  return { records, fileTypes };
}

function shouldPatchRecord(record, targetPatch) {
  return !targetPatch || targetPatch.records.has(record);
}

function shouldPatchFileType(fileType, targetPatch) {
  return !targetPatch || targetPatch.fileTypes.has(fileType);
}

function skipBlock(lines, startIndex, valueLineCount) {
  let index = startIndex;
  let valueRows = 0;

  while (index + 1 < lines.length && valueRows < valueLineCount) {
    index++;
    if (tokenizeEPolyScatLine(lines[index]).length > 0) {
      valueRows++;
    }
  }

  return index;
}

function appendMissingRecords(lines, seenRecords, replacements, targetPatch) {
  SINGLE_LINE_RECORD_ORDER.forEach(record => {
    if (seenRecords.has(record) || !replacements[record] || !shouldPatchRecord(record, targetPatch)) {
      return;
    }

    lines.push(replacements[record]);
  });
}

function appendMissingFileNames(lines, seenFileTypes, fileNameReplacements, targetPatch) {
  Object.keys(fileNameReplacements).forEach(fileType => {
    if (seenFileTypes.has(fileType) || !shouldPatchFileType(fileType, targetPatch)) {
      return;
    }

    lines.push(fileNameReplacements[fileType]);
  });
}

export function patchEPolyScatInputScript(contents, values, outputDefinitions = [], options = {}) {
  if (!String(contents || "").trim()) {
    return buildEPolyScatInputScript(values, outputDefinitions);
  }

  const targetPatch = getTargetPatch(options, outputDefinitions);
  const replacements = buildSingleLineReplacements(values);
  const fileNameReplacements = buildFileNameReplacements(values, outputDefinitions);
  const lines = String(contents).split(/\r?\n/);
  const patchedLines = [];
  const seenRecords = new Set();
  const seenFileTypes = new Set();
  let titleSeen = false;

  for (let index = 0; index < lines.length; index++) {
    const line = lines[index];
    const trimmedLine = line.trim();

    if (trimmedLine.startsWith("#") && !titleSeen) {
      titleSeen = true;
      if (hasValue(values.title) && shouldPatchRecord("Title", targetPatch)) {
        patchedLines.push(`# ${values.title}`);
      } else {
        patchedLines.push(line);
      }
      continue;
    }

    const tokens = tokenizeEPolyScatLine(line);
    if (tokens.length === 0) {
      patchedLines.push(line);
      continue;
    }

    const record = normalizeEPolyScatRecord(tokens[0]);
    if (record === "EngForm" && shouldPatchRecord(record, targetPatch)) {
      patchedLines.push(buildEngFormBlock(values));
      seenRecords.add(record);
      index = skipBlock(lines, index, 2 + (Number(values.engFormTerms) || 0));
      continue;
    }

    if (record === "AsyPol" && shouldPatchRecord(record, targetPatch)) {
      patchedLines.push(buildAsyPolBlock(values));
      seenRecords.add(record);
      index = skipBlock(lines, index, 7);
      continue;
    }

    if (record === "FileName") {
      const fileType = tokens[1];
      if (fileNameReplacements[fileType] && shouldPatchFileType(fileType, targetPatch)) {
        patchedLines.push(fileNameReplacements[fileType]);
        seenFileTypes.add(fileType);
      } else {
        patchedLines.push(line);
      }
      continue;
    }

    if (replacements[record] && shouldPatchRecord(record, targetPatch)) {
      patchedLines.push(replacements[record]);
      seenRecords.add(record);
      continue;
    }

    patchedLines.push(line);
  }

  if (!titleSeen && hasValue(values.title) && shouldPatchRecord("Title", targetPatch)) {
    patchedLines.unshift(`# ${values.title}`);
  }

  appendMissingRecords(patchedLines, seenRecords, replacements, targetPatch);
  appendMissingFileNames(patchedLines, seenFileTypes, fileNameReplacements, targetPatch);

  return patchedLines.join("\n").replace(/\n{3,}/g, "\n\n");
}
