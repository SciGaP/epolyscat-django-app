const DATA_RECORD_LABELS = [
  "AsyPol",
  "CnvOrbSel",
  "DPotEng",
  "EMax",
  "EngForm",
  "EpsAsym",
  "FegeEng",
  "FixCharge",
  "GrnType",
  "IPot",
  "InitSpinDeg",
  "InitSym",
  "IterMax",
  "JMax",
  "JPPMax",
  "LFTimeDelayAngles",
  "LMax",
  "LMaxA",
  "LMaxI",
  "LMaxK",
  "Label",
  "MFTimeDelayAngles",
  "NECenter",
  "OrbOcc",
  "OrbOccInit",
  "PCutRd",
  "PlaneWvCharge",
  "PosFile",
  "PosFitL",
  "PosGridTol",
  "PosPlot",
  "PrintBlm",
  "ResSearchEng",
  "RotConstants",
  "ScatContSym",
  "ScatEng",
  "ScatEngN",
  "ScatSym",
  "SpinDeg",
  "TargSpinDeg",
  "TargSym",
  "TestOut",
  "VCorr",
  "VibAveNInp",
  "ViewOrbGrid",
  // FileName is documented with the commands, but it mutates ordered output state.
  "FileName",
];

const MULTILINE_DATA_RECORDS = new Set([
  "AsyPol",
  "CnvOrbSel",
  "EngForm",
  "LFTimeDelayAngles",
  "MFTimeDelayAngles",
  "OrbOcc",
  "OrbOccInit",
  "RotConstants",
  "VibAveNInp",
  "ViewOrbGrid",
]);

const COMMAND_LABELS = [
  "CalcInt",
  "CatFiles",
  "Convert",
  "DipoleOp",
  "DumpBasis",
  "DumpDataRecords",
  "DumpIdy",
  "DumpIdyAll",
  "DumpMesa",
  "DumpOrb",
  "DumpRecords",
  "EDCS",
  "EngToler",
  "Exit",
  "ExpOrb",
  "FindQ",
  "GenFormPhIon",
  "GenFormScat",
  "GeomNormMode",
  "GetBlms",
  "GetCro",
  "GetDPot",
  "GetPot",
  "LFTimeDelay",
  "MatrixElementsCombine",
  "MatrixElementsCollect",
  "MFDCS",
  "MFTimeDelay",
  "OrientCro",
  "OrientNCro",
  "PhIon",
  "PhIonN",
  "PhIonPlaneWv",
  "PrintFlag",
  "ReadBlms",
  "ResSearch",
  "ResWvFun",
  "RmDataRecord",
  "RotOrientAsym",
  "SaveBlms",
  "SaveTMats",
  "Scat",
  "ScatN",
  "ScatPos",
  "SchmidtOrth",
  "SymNormMode",
  "TotalCrossSection",
  "VibAve",
  "VibAveN",
  "VibAveNI",
  "ViewOrb",
];

const MULTILINE_COMMANDS = new Set([
  "GetCro",
  "LFTimeDelay",
  "MFTimeDelay",
]);

const DATA_RECORD_UI = {
  LMax: { fieldKey: "lMax", group: "grid-expansion", valueType: "number" },
  LMaxI: { fieldKey: "lMaxI", group: "grid-expansion", valueType: "number" },
  EMax: { fieldKey: "eMax", group: "grid-expansion", valueType: "number" },
  FegeEng: { fieldKey: "fegeEng", group: "energies", valueType: "number" },
  ScatEng: { fieldKey: "scatEng", group: "energies", valueType: "list" },
  LMaxK: { fieldKey: "lMaxK", group: "energies", valueType: "number" },
  VCorr: { fieldKey: "vCorr", group: "potentials", valueType: "string" },
  ScatContSym: { fieldKey: "scatContSym", group: "state-definitions", valueType: "string" },
  ScatSym: { fieldKey: "scatSym", group: "state-definitions", valueType: "string" },
  InitSym: { fieldKey: "initSym", group: "state-definitions", valueType: "string" },
  InitSpinDeg: { fieldKey: "initSpinDeg", group: "state-definitions", valueType: "number" },
  OrbOccInit: { fieldKey: "orbOccInit", group: "state-definitions", valueType: "list" },
  OrbOcc: { fieldKey: "orbOcc", group: "state-definitions", valueType: "list" },
  SpinDeg: { fieldKey: "spinDeg", group: "state-definitions", valueType: "number" },
  TargSym: { fieldKey: "targSym", group: "state-definitions", valueType: "string" },
  TargSpinDeg: { fieldKey: "targSpinDeg", group: "state-definitions", valueType: "number" },
  IPot: { fieldKey: "iPot", group: "energies", valueType: "number" },
  EngForm: { group: "potentials", continuation: "conditional" },
  AsyPol: { group: "potentials", continuation: "conditional" },
  FileName: { group: "outputs", repeatable: true },
};

export const EPOLYSCAT_INPUT_SCHEMA = Object.freeze({
  dataRecords: DATA_RECORD_LABELS.map(label => Object.freeze({
    label,
    aliases: [label.toLowerCase()],
    repeatable: true,
    continuation: MULTILINE_DATA_RECORDS.has(label) ? "until-keyword" : "none",
    ...(DATA_RECORD_UI[label] || {}),
    provenance: "ePolyScat manual and official sample jobs",
  })),
  commands: COMMAND_LABELS.map(label => Object.freeze({
    label,
    aliases: [label.toLowerCase()],
    repeatable: true,
    continuation: MULTILINE_COMMANDS.has(label) ? "until-keyword" : "none",
    provenance: "ePolyScat manual command index",
  })),
});

export const EPOLYSCAT_RECOMMENDED_VALUES = Object.freeze({
  title: "electron scattering from CH4 in A1 symmetry",
  calculationKind: "scattering",
  lMax: 15,
  lMaxI: 120,
  eMax: 50.0,
  lMaxK: 4,
  initSym: "SG",
  initSpinDeg: 1,
  orbOccInit: "2 2 2 2 2 4",
  orbOcc: "2 2 2 2 1 4",
  spinDeg: 1,
  targSym: "SG",
  targSpinDeg: 2,
  iPot: 15.581,
  engFormCharge: 0,
  engFormType: 1,
  engFormTerms: 3,
  vCorr: "PZ",
  asyPolSwitchD: 0.15,
  asyPolTerms: 1,
  asyPolCenter: 1,
  asyPolValue: 17.50,
  fegeEng: 13.0,
  scatContSym: "A1",
  scatSym: "SU",
  scatEng: "0.0001 0.01 0.5",
  convertSource: "$pt/input.molden",
  convertFormat: "molden",
  matrixElementsFile: "matrix-elements.idy",
  plotDataFile: "cross-sections.dat",
  aWaveFunFile: "",
  debugOutStem: "",
  dumpOrbBasisFile: "",
  dumpOutFile: "",
  mfdcsGeomFile: "",
  mfdcsFile: "",
  mfdcsFullFile: "",
  orientAsymDataFile: "",
  orientAsymEigFile: "",
  orientAsymGeomFile: "",
  orientDataFile: "",
  orientGeomFile: "",
  plotData2DFile: "",
  sWaveFunFile: "",
  vibAveIdyFile: "",
  viewOrbFile: "",
  viewOrbFluxFile: "",
  viewOrbGeomFile: "",
});

export const EPOLYSCAT_DATA_ENTRY_SECTIONS = Object.freeze([
  {
    id: "grid-expansion",
    label: "Grid / Expansion",
    recordGroup: "Basic grid records",
    type: "fields",
    fields: [
      { key: "title", label: "Title", record: "#", inputType: "text" },
      {
        key: "calculationKind",
        label: "Calculation kind",
        record: "Flow",
        inputType: "select",
        options: [
          { value: "scattering", text: "Electron scattering" },
          { value: "photoionization", text: "Photoionization" },
        ],
      },
      { key: "lMax", label: "LMax", record: "LMax", inputType: "number" },
      { key: "lMaxI", label: "LMaxI", record: "LMaxI", inputType: "number" },
      { key: "eMax", label: "EMax", record: "EMax", inputType: "number", step: "0.1" },
      { key: "lMaxK", label: "LMaxK", record: "LMaxK", inputType: "number" },
      { key: "convertSource", label: "Convert source", record: "Convert", inputType: "text" },
      { key: "convertFormat", label: "Convert format", record: "Convert", inputType: "text" },
    ],
  },
  {
    id: "state-definitions",
    label: "State Definitions",
    recordGroup: "Initial and target state records",
    type: "fields",
    fields: [
      { key: "initSym", label: "InitSym", record: "InitSym", inputType: "text" },
      { key: "initSpinDeg", label: "InitSpinDeg", record: "InitSpinDeg", inputType: "number" },
      { key: "orbOccInit", label: "OrbOccInit", record: "OrbOccInit", inputType: "text" },
      { key: "orbOcc", label: "OrbOcc", record: "OrbOcc", inputType: "text" },
      { key: "spinDeg", label: "SpinDeg", record: "SpinDeg", inputType: "number" },
      { key: "targSym", label: "TargSym", record: "TargSym", inputType: "text" },
      { key: "targSpinDeg", label: "TargSpinDeg", record: "TargSpinDeg", inputType: "number" },
      { key: "iPot", label: "IPot", record: "IPot", inputType: "number", step: "0.001" },
    ],
  },
  {
    id: "potentials",
    label: "Potentials",
    recordGroup: "EngForm, VCorr, AsyPol, FEGE",
    type: "fields",
    fields: [
      { key: "engFormCharge", label: "EngForm charge", record: "EngForm", inputType: "number" },
      { key: "engFormType", label: "EngForm type", record: "EngForm", inputType: "number" },
      { key: "engFormTerms", label: "EngForm terms", record: "EngForm", inputType: "number" },
      { key: "vCorr", label: "VCorr", record: "VCorr", inputType: "text" },
      { key: "asyPolSwitchD", label: "SwitchD", record: "AsyPol", inputType: "number", step: "0.01" },
      { key: "asyPolTerms", label: "nterm", record: "AsyPol", inputType: "number" },
      { key: "asyPolCenter", label: "Center", record: "AsyPol", inputType: "number" },
      { key: "asyPolValue", label: "Polarizability", record: "AsyPol", inputType: "number", step: "0.01" },
      { key: "fegeEng", label: "FegeEng", record: "FegeEng", inputType: "number", step: "0.1" },
    ],
  },
  {
    id: "energies-partial-waves",
    label: "Energies / Partial Waves",
    recordGroup: "Scattering symmetry and energy records",
    type: "fields",
    fields: [
      { key: "scatContSym", label: "ScatContSym", record: "ScatContSym", inputType: "text" },
      { key: "scatSym", label: "ScatSym", record: "ScatSym", inputType: "text" },
      { key: "scatEng", label: "ScatEng", record: "ScatEng", inputType: "text" },
    ],
  },
  {
    id: "outputs",
    label: "Outputs",
    recordGroup: "FileName records",
    type: "outputs",
  },
]);

export const EPOLYSCAT_OUTPUT_DEFINITIONS = Object.freeze([
  {
    fileType: "AWaveFun",
    valueKey: "aWaveFunFile",
    extension: ".dat",
    description: "Adiabatic resonance wave functions",
    disposition: "REWIND",
  },
  {
    fileType: "DebugOut",
    valueKey: "debugOutStem",
    extension: "",
    description: "Per-process debugging output stem",
    disposition: "REWIND",
  },
  {
    fileType: "DumpOrbBasis",
    valueKey: "dumpOrbBasisFile",
    extension: ".dat",
    description: "DumpOrb and DumpBasis output",
    disposition: "REWIND",
  },
  {
    fileType: "DumpOut",
    valueKey: "dumpOutFile",
    extension: ".dat",
    description: "DumpIdy output",
    disposition: "REWIND",
  },
  {
    fileType: "MatrixElements",
    valueKey: "matrixElementsFile",
    extension: ".idy",
    description: "Dynamical coefficients",
    disposition: "REWIND",
  },
  {
    fileType: "MFDCSGeom",
    valueKey: "mfdcsGeomFile",
    extension: ".dat",
    description: "Single-orientation MFDCS geometry",
    disposition: "REWIND",
  },
  {
    fileType: "MFDCS",
    valueKey: "mfdcsFile",
    extension: ".dat",
    description: "Two-dimensional differential cross sections",
    disposition: "REWIND",
  },
  {
    fileType: "MFDCSFull",
    valueKey: "mfdcsFullFile",
    extension: ".dat",
    description: "Full-list differential cross sections",
    disposition: "REWIND",
  },
  {
    fileType: "OrientAsymData",
    valueKey: "orientAsymDataFile",
    extension: ".dat",
    description: "Rotated asymmetric-top molecular-frame data",
    disposition: "REWIND",
  },
  {
    fileType: "OrientAsymEig",
    valueKey: "orientAsymEigFile",
    extension: ".dat",
    description: "Asymmetric-top rotational eigenvalues",
    disposition: "REWIND",
  },
  {
    fileType: "OrientAsymGeom",
    valueKey: "orientAsymGeomFile",
    extension: ".dat",
    description: "Asymmetric-top geometry data",
    disposition: "REWIND",
  },
  {
    fileType: "OrientData",
    valueKey: "orientDataFile",
    extension: ".dat",
    description: "Oriented molecule analysis",
    disposition: "REWIND",
  },
  {
    fileType: "OrientGeom",
    valueKey: "orientGeomFile",
    extension: ".dat",
    description: "Geometry for oriented cross sections",
    disposition: "REWIND",
  },
  {
    fileType: "PlotData",
    valueKey: "plotDataFile",
    extension: ".dat",
    description: "One-variable plot data",
    disposition: "REWIND",
  },
  {
    fileType: "PlotData2D",
    valueKey: "plotData2DFile",
    extension: ".dat",
    description: "Two-variable plot data",
    disposition: "REWIND",
  },
  {
    fileType: "SWaveFun",
    valueKey: "sWaveFunFile",
    extension: ".dat",
    description: "Spherical resonance wave functions",
    disposition: "REWIND",
  },
  {
    fileType: "VibAveIdy",
    valueKey: "vibAveIdyFile",
    extension: ".idy",
    description: "Vibrationally averaged dipole matrix elements",
    disposition: "REWIND",
  },
  {
    fileType: "ViewOrb",
    valueKey: "viewOrbFile",
    extension: ".dat",
    description: "Orbital or potential values",
    disposition: "REWIND",
  },
  {
    fileType: "ViewOrbFlux",
    valueKey: "viewOrbFluxFile",
    extension: ".dat",
    description: "Orbital flux and divergence",
    disposition: "REWIND",
  },
  {
    fileType: "ViewOrbGeom",
    valueKey: "viewOrbGeomFile",
    extension: ".dat",
    description: "Geometry and grid values",
    disposition: "REWIND",
  },
]);

const DATA_RECORD_BY_TOKEN = EPOLYSCAT_INPUT_SCHEMA.dataRecords.reduce((records, record) => {
  [record.label, ...(record.aliases || [])].forEach(alias => {
    records[String(alias).toLowerCase()] = record;
  });
  return records;
}, {});

const COMMAND_BY_TOKEN = EPOLYSCAT_INPUT_SCHEMA.commands.reduce((commands, command) => {
  [command.label, ...(command.aliases || [])].forEach(alias => {
    commands[String(alias).toLowerCase()] = command;
  });
  return commands;
}, {});

function splitSourceLines(contents) {
  const source = String(contents == null ? "" : contents);
  const lines = [];
  const pattern = /([^\r\n]*)(\r\n|\n|\r|$)/g;
  let match;
  let offset = 0;
  let lineNumber = 1;

  while ((match = pattern.exec(source)) !== null) {
    if (match[0] === "") {
      break;
    }
    lines.push({
      raw: match[1],
      eol: match[2],
      start: offset,
      end: offset + match[0].length,
      line: lineNumber,
    });
    offset += match[0].length;
    lineNumber += 1;
  }
  return lines;
}

function sourceLinesText(sourceLines) {
  return sourceLines.map(line => `${line.raw}${line.eol}`).join("");
}

function makeSourceSpan(sourceLines) {
  const firstLine = sourceLines[0];
  const lastLine = sourceLines[sourceLines.length - 1];
  return {
    start: firstLine.start,
    end: lastLine.end,
    startLine: firstLine.line,
    endLine: lastLine.line,
  };
}

function firstLineKeyword(line) {
  const tokens = tokenizeEPolyScatLine(line.raw);
  return tokens.length ? tokens[0] : "";
}

function engFormBlockIsComplete(sourceLines) {
  const structure = engFormSourceStructure(sourceLines);
  if (!structure.formulaRow || structure.formulaRow.tokens.length < 2) {
    return false;
  }
  if (structure.formulaType === 0) {
    return true;
  }
  if (![1, 2].includes(structure.formulaType)) {
    return false;
  }
  return Boolean(structure.countRow)
      && Number.isInteger(structure.termCount)
      && structure.termRows.length >= structure.termCount;
}

function isContinuationLine(line, dataRecord, sourceLines) {
  if (dataRecord.label === "EngForm" && engFormBlockIsComplete(sourceLines)) {
    return false;
  }
  if (line.raw.trim().startsWith("#")) {
    return true;
  }
  if (!line.raw.trim()) {
    return true;
  }
  const tokens = tokenizeEPolyScatLine(line.raw);
  if (tokens.length === 0) {
    return false;
  }
  const firstToken = String(tokens[0]);
  if (DATA_RECORD_BY_TOKEN[firstToken.toLowerCase()] || COMMAND_BY_TOKEN[firstToken.toLowerCase()]) {
    return false;
  }
  return true;
}

function isCommandContinuationLine(line) {
  const trimmed = line.raw.trim();
  if (!trimmed || trimmed.startsWith("#")) {
    return false;
  }
  const keyword = firstLineKeyword(line).toLowerCase();
  return !DATA_RECORD_BY_TOKEN[keyword] && !COMMAND_BY_TOKEN[keyword];
}

function buildDocumentNode(type, sourceLines, extra = {}) {
  return {
    type,
    raw: sourceLinesText(sourceLines),
    sourceLines,
    sourceSpan: makeSourceSpan(sourceLines),
    ...extra,
  };
}

export function parseEPolyScatDocument(contents) {
  const source = String(contents == null ? "" : contents);
  const lines = splitSourceLines(source);
  const nodes = [];

  for (let index = 0; index < lines.length; index++) {
    const line = lines[index];
    const trimmed = line.raw.trim();
    if (!trimmed) {
      nodes.push(buildDocumentNode("BlankLine", [line]));
      continue;
    }
    if (trimmed.startsWith("#")) {
      nodes.push(buildDocumentNode("Comment", [line], {
        text: trimmed.replace(/^#+\s*/, ""),
      }));
      continue;
    }

    const keyword = firstLineKeyword(line);
    const dataRecord = DATA_RECORD_BY_TOKEN[String(keyword).toLowerCase()];
    const command = COMMAND_BY_TOKEN[String(keyword).toLowerCase()];
    if (dataRecord) {
      const sourceLines = [line];
      const firstLineTokens = tokenizeEPolyScatLine(line.raw);
      if (dataRecord.continuation !== "none" || firstLineTokens.length === 1) {
        while (
          index + 1 < lines.length
          && isContinuationLine(lines[index + 1], dataRecord, sourceLines)
        ) {
          sourceLines.push(lines[index + 1]);
          index += 1;
        }
      }
      const tokens = firstLineTokens;
      nodes.push(buildDocumentNode("DataRecord", sourceLines, {
        label: dataRecord.label,
        arguments: tokens.slice(1),
        continuationRows: sourceLines.slice(1).map(sourceLine => ({
          raw: sourceLine.raw,
          tokens: tokenizeEPolyScatLine(sourceLine.raw),
        })),
        schema: dataRecord,
      }));
      continue;
    }
    if (command) {
      const sourceLines = [line];
      if (command.continuation !== "none") {
        while (
          index + 1 < lines.length
          && isCommandContinuationLine(lines[index + 1])
        ) {
          sourceLines.push(lines[index + 1]);
          index += 1;
        }
      }
      const tokens = tokenizeEPolyScatLine(line.raw);
      nodes.push(buildDocumentNode("Command", sourceLines, {
        label: command.label,
        arguments: tokens.slice(1),
        continuationRows: sourceLines.slice(1).map(sourceLine => ({
          raw: sourceLine.raw,
          tokens: tokenizeEPolyScatLine(sourceLine.raw),
        })),
        schema: command,
      }));
      continue;
    }
    nodes.push(buildDocumentNode("Unknown", [line], { keyword }));
  }

  const newlineLine = lines.find(line => line.eol);
  return {
    type: "Document",
    source,
    newline: newlineLine ? newlineLine.eol : "\n",
    nodes,
  };
}

export function serializeEPolyScatDocument(document) {
  if (!document || !Array.isArray(document.nodes)) {
    return "";
  }
  return document.nodes.map(node => sourceLinesText(node.sourceLines || [])).join("");
}

const FIELD_RECORDS = {
  title: ["Title"],
  lMax: ["LMax"],
  lMaxI: ["LMaxI"],
  eMax: ["EMax"],
  fegeEng: ["FegeEng"],
  scatEng: ["ScatEng"],
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

  const token = String(record || "").toLowerCase();
  const schemaEntry = DATA_RECORD_BY_TOKEN[token] || COMMAND_BY_TOKEN[token];
  return recordMap[token] || (schemaEntry ? schemaEntry.label : record);
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

function nodeValueRows(node) {
  return (node.continuationRows || [])
      .map(row => row.tokens || tokenizeEPolyScatLine(row.raw || ""))
      .filter(tokens => tokens.length > 0);
}

function semanticSourceRows(sourceLines, startIndex = 0) {
  return sourceLines
      .map((sourceLine, lineIndex) => ({
        lineIndex,
        sourceLine,
        tokens: tokenizeEPolyScatLine(sourceLine.raw),
      }))
      .filter(row => row.lineIndex >= startIndex && row.tokens.length > 0);
}

function engFormSourceStructure(sourceLines) {
  const headerTokens = tokenizeEPolyScatLine(sourceLines[0] ? sourceLines[0].raw : "");
  const inlineFormula = headerTokens.length >= 3;
  const continuationRows = semanticSourceRows(sourceLines, 1);
  const formulaRow = inlineFormula
      ? {
        lineIndex: 0,
        sourceLine: sourceLines[0],
        tokens: headerTokens.slice(1),
        inline: true,
      }
      : continuationRows[0] || null;
  const formulaType = formulaRow && formulaRow.tokens.length >= 2
      ? Number(formulaRow.tokens[1])
      : NaN;
  const countRowIndex = inlineFormula ? 0 : 1;
  const countRow = [1, 2].includes(formulaType)
      ? continuationRows[countRowIndex] || null
      : null;
  const termCount = countRow ? Number(countRow.tokens[0]) : 0;
  const termRows = countRow && Number.isInteger(termCount) && termCount >= 0
      ? continuationRows.slice(countRowIndex + 1, countRowIndex + 1 + termCount)
      : [];

  return {
    inlineFormula,
    formulaRow,
    formulaType,
    countRow,
    termCount,
    termRows,
  };
}

function asyPolStructure(node) {
  const continuationRows = semanticSourceRows(node.sourceLines, 1);
  const rows = (node.arguments || []).length
      ? [{
        lineIndex: 0,
        sourceLine: node.sourceLines[0],
        tokens: node.arguments,
        inline: true,
      }, ...continuationRows]
      : continuationRows;
  const switchRow = rows[0] || null;
  const termCountRow = rows[1] || null;
  const centerRow = rows[2] || null;
  let cursor = 3;
  const center = centerRow ? Number(centerRow.tokens[0]) : NaN;
  const centerCoordinatesRow = center === 0 ? rows[cursor++] || null : null;
  const typeRow = rows[cursor++] || null;
  const valueRow = rows[cursor] || null;

  return {
    switchRow,
    termCountRow,
    centerRow,
    centerCoordinatesRow,
    typeRow,
    valueRow,
  };
}

function parseEngFormNode(node) {
  const structure = engFormSourceStructure(node.sourceLines);
  const values = {};

  if (structure.formulaRow) {
    setParsedValue(values, "engFormCharge", parseNumericToken(structure.formulaRow.tokens[0]));
    setParsedValue(values, "engFormType", parseNumericToken(structure.formulaRow.tokens[1]));
  }
  if (structure.countRow && Number(values.engFormType) !== 0) {
    setParsedValue(values, "engFormTerms", parseNumericToken(structure.countRow.tokens[0]));
  }
  return values;
}

function parseAsyPolNode(node) {
  const structure = asyPolStructure(node);
  const values = {};

  if (structure.switchRow) {
    setParsedValue(values, "asyPolSwitchD", parseNumericToken(structure.switchRow.tokens[0]));
  }
  if (structure.termCountRow) {
    setParsedValue(values, "asyPolTerms", parseNumericToken(structure.termCountRow.tokens[0]));
  }
  if (structure.centerRow) {
    setParsedValue(values, "asyPolCenter", parseNumericToken(structure.centerRow.tokens[0]));
  }
  if (structure.valueRow) {
    setParsedValue(values, "asyPolValue", parseNumericToken(structure.valueRow.tokens[0]));
  }
  return values;
}

function recordTokens(node) {
  return [node.label, ...(node.arguments || [])];
}

export function parseEPolyScatInputScript(contents, outputDefinitions = []) {
  const normalizedContents = normalizeEPolyScatInputContents(contents);
  const document = parseEPolyScatDocument(normalizedContents);
  const parsedValues = {};

  document.nodes.forEach(node => {
    if (node.type === "Comment" && !parsedValues.title && node.text) {
      parsedValues.title = node.text.trim();
      return;
    }
    if (node.type !== "DataRecord" && node.type !== "Command") {
      return;
    }

    const record = node.label;
    const tokens = recordTokens(node);
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
      const values = tokens.slice(1).length
          ? tokens.slice(1)
          : nodeValueRows(node).reduce((items, row) => items.concat(row), []);
      parsedValues.orbOccInit = values.join(" ");
      parsedValues.calculationKind = "photoionization";
    } else if (record === "OrbOcc") {
      const values = tokens.slice(1).length
          ? tokens.slice(1)
          : nodeValueRows(node).reduce((items, row) => items.concat(row), []);
      parsedValues.orbOcc = values.join(" ");
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
      Object.assign(parsedValues, parseEngFormNode(node));
    } else if (record === "AsyPol") {
      Object.assign(parsedValues, parseAsyPolNode(node));
    } else if (["GenFormPhIon", "DipoleOp", "PhIon", "PhIonN", "GetCro"].includes(record)) {
      parsedValues.calculationKind = "photoionization";
    }
  });

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
  const formulaType = Number(values.engFormType) || 0;
  const header = [
    "EngForm",
    ` ${values.engFormCharge} ${formulaType}`,
  ];
  if (formulaType === 0) {
    return header.join("\n");
  }

  const termCount = Number(values.engFormTerms) || 0;
  const defaultTerm = formulaType === 2 ? " 2.0 -1.0 0" : " 2.0 -1.0";
  const terms = Array.isArray(values.engFormRows)
      ? values.engFormRows.map(row => ` ${Array.isArray(row) ? row.join(" ") : row}`)
      : Array.from({ length: termCount }, () => defaultTerm);
  return [
    ...header,
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

  const generated = lines.join("\n").replace(/\n{3,}/g, "\n\n");
  return serializeEPolyScatDocument(parseEPolyScatDocument(generated));
}

function buildSingleLineReplacements(values) {
  const replacements = {};

  if (hasValue(values.lMax)) replacements.LMax = `LMax ${values.lMax}`;
  if (hasValue(values.lMaxI)) replacements.LMaxI = `LMaxI ${values.lMaxI}`;
  if (hasValue(values.eMax)) replacements.EMax = `EMax ${values.eMax}`;
  if (hasValue(values.fegeEng)) replacements.FegeEng = `FegeEng ${values.fegeEng}`;
  if (hasValue(values.scatEng)) {
    replacements.ScatEng = `ScatEng ${values.scatEng}`;
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

  return { records, fileTypes, changedKeys: new Set(changedKeys) };
}

function shouldPatchRecord(record, targetPatch) {
  return !targetPatch || targetPatch.records.has(record);
}

function shouldPatchFileType(fileType, targetPatch) {
  return !targetPatch || targetPatch.fileTypes.has(fileType);
}

function shouldPatchField(fieldKey, targetPatch) {
  return !targetPatch || targetPatch.changedKeys.has(fieldKey);
}

function refreshDocumentNode(node) {
  node.raw = sourceLinesText(node.sourceLines);
  if (node.type === "DataRecord" || node.type === "Command") {
    const tokens = tokenizeEPolyScatLine(node.sourceLines[0].raw);
    node.arguments = tokens.slice(1);
    node.continuationRows = node.sourceLines.slice(1).map(sourceLine => ({
      raw: sourceLine.raw,
      tokens: tokenizeEPolyScatLine(sourceLine.raw),
    }));
  }
}

function replaceSourceLineCode(sourceLine, replacement) {
  const oldCode = stripEPolyScatComment(sourceLine.raw);
  const comment = sourceLine.raw.slice(oldCode.length);
  const indentation = (sourceLine.raw.match(/^\s*/) || [""])[0];
  const trailingWhitespace = (oldCode.match(/\s*$/) || [""])[0];
  sourceLine.raw = `${indentation}${String(replacement).trim()}${comment ? trailingWhitespace + comment : ""}`;
}

function replaceNodeCode(node, replacement) {
  replaceSourceLineCode(node.sourceLines[0], replacement);
  refreshDocumentNode(node);
}

function lastMatchingNode(document, predicate) {
  for (let index = document.nodes.length - 1; index >= 0; index--) {
    if (predicate(document.nodes[index])) {
      return document.nodes[index];
    }
  }
  return null;
}

function nodeMatchesLabel(node, label) {
  return (node.type === "DataRecord" || node.type === "Command") && node.label === label;
}

function insertSourceBeforeCommands(document, source) {
  const newline = document.newline || "\n";
  const insertionIndex = document.nodes.findIndex(node => node.type === "Command");
  const targetIndex = insertionIndex >= 0 ? insertionIndex : document.nodes.length;
  const previousNode = document.nodes[targetIndex - 1];

  if (previousNode && previousNode.sourceLines.length) {
    const previousLine = previousNode.sourceLines[previousNode.sourceLines.length - 1];
    if (!previousLine.eol) {
      previousLine.eol = newline;
      refreshDocumentNode(previousNode);
    }
  }

  const inserted = parseEPolyScatDocument(`${String(source).replace(/[\r\n]+$/, "")}${newline}`);
  document.nodes.splice(targetIndex, 0, ...inserted.nodes);
}

function removeNodeSourceRows(node, rows) {
  const sourceLines = new Set(rows.filter(Boolean).map(row => row.sourceLine));
  node.sourceLines = node.sourceLines.filter(sourceLine => !sourceLines.has(sourceLine));
}

function insertNodeSourceLineAfter(node, anchorSourceLine, raw) {
  const newlineSourceLine = node.sourceLines.find(sourceLine => sourceLine.eol);
  const newline = newlineSourceLine ? newlineSourceLine.eol : "\n";
  const anchorIndex = node.sourceLines.indexOf(anchorSourceLine);
  if (anchorIndex < 0) {
    return null;
  }
  if (!anchorSourceLine.eol) {
    anchorSourceLine.eol = newline;
  }
  const sourceLine = splitSourceLines(`${raw}${newline}`)[0];
  node.sourceLines.splice(anchorIndex + 1, 0, sourceLine);
  return sourceLine;
}

function replaceSemanticRowValue(node, row, value) {
  if (!row) {
    return;
  }
  if (row.inline) {
    const headerTokens = tokenizeEPolyScatLine(node.sourceLines[0].raw);
    headerTokens[1] = String(value);
    replaceSourceLineCode(node.sourceLines[0], headerTokens.join(" "));
    return;
  }
  replaceSourceLineCode(row.sourceLine, value);
}

function patchEngFormNode(node, values, targetPatch) {
  if (!node || !shouldPatchRecord("EngForm", targetPatch)) {
    return false;
  }
  const originalStructure = engFormSourceStructure(node.sourceLines);
  const currentValues = parseEngFormNode(node);
  const patchCharge = shouldPatchField("engFormCharge", targetPatch);
  const patchType = shouldPatchField("engFormType", targetPatch);
  const patchTerms = shouldPatchField("engFormTerms", targetPatch);
  const charge = patchCharge && hasValue(values.engFormCharge)
      ? values.engFormCharge
      : currentValues.engFormCharge;
  const formulaType = patchType && hasValue(values.engFormType)
      ? values.engFormType
      : currentValues.engFormType;
  const hasFormulaUpdate = (patchCharge && hasValue(values.engFormCharge))
      || (patchType && hasValue(values.engFormType));
  const updateTermStructure = patchTerms || (patchType && hasValue(values.engFormType));

  if (!originalStructure.formulaRow || !hasValue(charge) || !hasValue(formulaType)) {
    return true;
  }

  if (hasFormulaUpdate) {
    if (originalStructure.inlineFormula) {
      replaceSourceLineCode(node.sourceLines[0], `EngForm ${charge} ${formulaType}`);
    } else {
      replaceSourceLineCode(originalStructure.formulaRow.sourceLine, `${charge} ${formulaType}`);
    }
  }

  if (!updateTermStructure) {
    refreshDocumentNode(node);
    return true;
  }

  const numericFormulaType = Number(formulaType);
  if (numericFormulaType === 0) {
    removeNodeSourceRows(node, [
      originalStructure.countRow,
      ...originalStructure.termRows,
    ]);
    refreshDocumentNode(node);
    return true;
  }
  if (![1, 2].includes(numericFormulaType)) {
    refreshDocumentNode(node);
    return true;
  }

  const desiredTermCount = hasValue(values.engFormTerms)
      ? Math.max(0, Number(values.engFormTerms) || 0)
      : Math.max(0, originalStructure.termCount || 0);
  let countSourceLine = originalStructure.countRow
      ? originalStructure.countRow.sourceLine
      : null;
  if (countSourceLine) {
    replaceSourceLineCode(countSourceLine, desiredTermCount);
  } else {
    countSourceLine = insertNodeSourceLineAfter(
        node,
        originalStructure.formulaRow.sourceLine,
        ` ${desiredTermCount}`
    );
  }

  const retainedTermRows = originalStructure.termRows.slice(0, desiredTermCount);
  removeNodeSourceRows(node, originalStructure.termRows.slice(desiredTermCount));
  retainedTermRows.forEach(row => {
    const tokens = [...row.tokens];
    if (numericFormulaType === 2 && tokens.length < 3) {
      tokens.push("0");
    } else if (numericFormulaType === 1 && tokens.length > 2) {
      tokens.splice(2);
    }
    replaceSourceLineCode(row.sourceLine, tokens.join(" "));
  });

  let anchorSourceLine = retainedTermRows.length
      ? retainedTermRows[retainedTermRows.length - 1].sourceLine
      : countSourceLine;
  const defaultTerm = numericFormulaType === 2 ? " 2.0 -1.0 0" : " 2.0 -1.0";
  for (let index = retainedTermRows.length; index < desiredTermCount; index++) {
    anchorSourceLine = insertNodeSourceLineAfter(node, anchorSourceLine, defaultTerm);
  }
  refreshDocumentNode(node);
  return true;
}

function patchAsyPolNode(node, values, targetPatch) {
  if (!node || !shouldPatchRecord("AsyPol", targetPatch)) {
    return false;
  }
  const structure = asyPolStructure(node);
  const updates = [
    ["asyPolSwitchD", structure.switchRow],
    ["asyPolTerms", structure.termCountRow],
    ["asyPolCenter", structure.centerRow],
    ["asyPolValue", structure.valueRow],
  ];

  updates.forEach(([key, row]) => {
    if (!shouldPatchField(key, targetPatch) || !hasValue(values[key]) || !row) {
      return;
    }
    replaceSemanticRowValue(node, row, values[key]);
  });
  refreshDocumentNode(node);
  return true;
}

function patchTitle(document, title, targetPatch) {
  if (!hasValue(title) || !shouldPatchRecord("Title", targetPatch)) {
    return;
  }
  const titleNode = document.nodes.find(node => node.type === "Comment" && node.text);
  if (titleNode) {
    replaceNodeCode(titleNode, `# ${title}`);
  } else {
    const newline = document.newline || "\n";
    const titleDocument = parseEPolyScatDocument(`# ${title}${newline}`);
    document.nodes.unshift(...titleDocument.nodes);
  }
}

function patchSingleLineRecords(document, replacements, targetPatch) {
  SINGLE_LINE_RECORD_ORDER.forEach(record => {
    const replacement = replacements[record];
    if (!replacement || !shouldPatchRecord(record, targetPatch)) {
      return;
    }
    const node = lastMatchingNode(document, item => nodeMatchesLabel(item, record));
    if (node) {
      replaceNodeCode(node, replacement);
    } else {
      insertSourceBeforeCommands(document, replacement);
    }
  });
}

function patchStructuredRecords(document, values, targetPatch) {
  const engFormNode = lastMatchingNode(document, node => nodeMatchesLabel(node, "EngForm"));
  if (!patchEngFormNode(engFormNode, values, targetPatch)
      && shouldPatchRecord("EngForm", targetPatch)
      && hasValue(values.engFormCharge)
      && hasValue(values.engFormType)) {
    insertSourceBeforeCommands(document, buildEngFormBlock(values));
  }

  const asyPolNode = lastMatchingNode(document, node => nodeMatchesLabel(node, "AsyPol"));
  if (!patchAsyPolNode(asyPolNode, values, targetPatch)
      && shouldPatchRecord("AsyPol", targetPatch)
      && hasValue(values.asyPolSwitchD)) {
    insertSourceBeforeCommands(document, buildAsyPolBlock(values));
  }
}

function patchFileNameRecords(document, replacements, targetPatch) {
  Object.keys(replacements).forEach(fileType => {
    if (!shouldPatchFileType(fileType, targetPatch)) {
      return;
    }
    const node = lastMatchingNode(document, item => (
      nodeMatchesLabel(item, "FileName") && item.arguments[0] === fileType
    ));
    if (node) {
      replaceNodeCode(node, replacements[fileType]);
    } else {
      insertSourceBeforeCommands(document, replacements[fileType]);
    }
  });
}

export function patchEPolyScatInputScript(contents, values, outputDefinitions = [], options = {}) {
  if (!String(contents || "").trim()) {
    return buildEPolyScatInputScript(values, outputDefinitions);
  }

  const targetPatch = getTargetPatch(options, outputDefinitions);
  const replacements = buildSingleLineReplacements(values);
  const fileNameReplacements = buildFileNameReplacements(values, outputDefinitions);
  const document = parseEPolyScatDocument(String(contents));

  patchTitle(document, values.title, targetPatch);
  patchSingleLineRecords(document, replacements, targetPatch);
  patchStructuredRecords(document, values, targetPatch);
  patchFileNameRecords(document, fileNameReplacements, targetPatch);

  return serializeEPolyScatDocument(document);
}

const LEGACY_PAGE_NAMES = {
  "grid-expansion": "Grid / Expansion",
  "state-definitions": "State Definitions",
  potentials: "Potentials",
  energies: "Energies / Partial Waves",
  outputs: "Outputs",
  additional: "Additional Records",
};

function legacyRecordValue(contents, record) {
  const parsedValues = parseEPolyScatInputScript(contents);
  if (record.fieldKey && hasValue(parsedValues[record.fieldKey])) {
    return parsedValues[record.fieldKey];
  }
  const document = parseEPolyScatDocument(contents);
  const node = lastMatchingNode(document, item => nodeMatchesLabel(item, record.label));
  if (!node) {
    return "";
  }
  if (node.arguments && node.arguments.length) {
    return node.arguments.join(" ");
  }
  return (node.continuationRows || []).map(row => row.tokens.join(" ")).join("; ");
}

async function setLegacyRecordValue({ getContents, setContents }, record, value) {
  const contents = await getContents();
  if (record.fieldKey) {
    const patched = patchEPolyScatInputScript(
        contents,
        { [record.fieldKey]: value },
        [],
        { changedKeys: [record.fieldKey] }
    );
    await setContents(patched);
    return;
  }

  const document = parseEPolyScatDocument(contents);
  const node = lastMatchingNode(document, item => nodeMatchesLabel(item, record.label));
  const replacement = `${record.label}${hasValue(value) ? ` ${value}` : ""}`;
  if (node) {
    replaceNodeCode(node, replacement);
  } else {
    insertSourceBeforeCommands(document, replacement);
  }
  await setContents(serializeEPolyScatDocument(document));
}

function legacyTableField(adapter, record) {
  return {
    name: record.label,
    type: record.valueType === "number" ? "number" : "text",
    step: record.valueType === "number" ? "any" : undefined,
    get: async () => legacyRecordValue(await adapter.getContents(), record),
    set: async value => setLegacyRecordValue(adapter, record, value),
    issues: () => [],
  };
}

function commandRows(adapter) {
  return async () => {
    const document = parseEPolyScatDocument(await adapter.getContents());
    const commands = document.nodes.filter(node => node.type === "Command");
    return commands.map((command, commandIndex) => ({
      delete: async () => {
        const currentDocument = parseEPolyScatDocument(await adapter.getContents());
        const currentCommands = currentDocument.nodes.filter(node => node.type === "Command");
        const currentCommand = currentCommands[commandIndex];
        if (currentCommand) {
          currentDocument.nodes = currentDocument.nodes.filter(node => node !== currentCommand);
          await adapter.setContents(serializeEPolyScatDocument(currentDocument));
        }
      },
      cells: [{
        get: async () => command.sourceLines[0].raw.trim(),
        set: async value => {
          const currentDocument = parseEPolyScatDocument(await adapter.getContents());
          const currentCommands = currentDocument.nodes.filter(node => node.type === "Command");
          const currentCommand = currentCommands[commandIndex];
          if (currentCommand) {
            replaceNodeCode(currentCommand, value);
            await adapter.setContents(serializeEPolyScatDocument(currentDocument));
          }
        },
      }],
    }));
  };
}

function appendCommand(adapter, command) {
  return async () => {
    const contents = await adapter.getContents();
    const document = parseEPolyScatDocument(contents);
    const newline = document.newline || "\n";
    const lastNode = document.nodes[document.nodes.length - 1];
    if (lastNode && lastNode.sourceLines.length) {
      const lastLine = lastNode.sourceLines[lastNode.sourceLines.length - 1];
      if (!lastLine.eol) {
        lastLine.eol = newline;
        refreshDocumentNode(lastNode);
      }
    }
    const commandDocument = parseEPolyScatDocument(command);
    document.nodes.push(...commandDocument.nodes);
    await adapter.setContents(serializeEPolyScatDocument(document));
  };
}

export function createLegacyEPolyScatTableObject(adapter) {
  if (!adapter || typeof adapter.getContents !== "function" || typeof adapter.setContents !== "function") {
    throw new TypeError("The legacy ePolyScat table adapter requires getContents and setContents.");
  }

  const pageRecords = EPOLYSCAT_INPUT_SCHEMA.dataRecords.reduce((pages, record) => {
    const group = record.group || "additional";
    pages[group] = pages[group] || [];
    pages[group].push(record);
    return pages;
  }, {});
  const pages = Object.keys(LEGACY_PAGE_NAMES)
      .filter(group => pageRecords[group] && pageRecords[group].length)
      .map(group => ({
        name: LEGACY_PAGE_NAMES[group],
        type: "table",
        data: pageRecords[group].map(record => legacyTableField(adapter, record)),
      }));

  pages.push({
    name: "Command Sequence",
    type: "list",
    columns: [{ name: "ePolyScat command", type: "text", issues: () => [] }],
    get: commandRows(adapter),
    addRow: appendCommand(adapter, "GetBlms"),
  });
  return { pages };
}
