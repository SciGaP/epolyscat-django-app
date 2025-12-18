import store from "./store";

export const descriptions = {
    "MODULE": "Individual Steps in the ePolyScat workflow",
        "Gaussian16": "Computes inital wavefunctions for the ePolyScat  Calculations",
        "OpenMolcas": "Computes inital wavefunctions for the ePolyScat Calculations",
        "ePolyScat": "Computes scattering cross sections using the ePolyScat program",
    "UTILITY": "Small programs ment for preparing input/output data",
        "MoldenMerge": "Merge Molden data files into one file",
        "CnvMath": "Convert data into Mathematica format files",
        "CnvMatLab": "Converts data into Matlab files",
        "CnvLinFull": "Computes the Phonon ion differential cross section ",
        "NRFPAD": "computes the final observable N photon RFPAD from the differential cross section for an oriented molecular calculation",
        "Cube2igor": "Converts the G16 Cube data files into IGOR plot files",
    "WORKFLOW": "Pre-constructed sequences of modules and utilities designed to carry out entire processes",
        "Data_Gen": "Data_Generation",
        "ePolyScat_Run": "ePolyScat_Run",
        "Analysis": "Analysis",
        "BOUND": "Does bound state calculations, produces bound_tab file",
        "OSCPOL": "OSCPOL",
        "STGF": "STGF",
        "SCATTERING": "produces H.DAT file",
    "Intermediate Files": "These are all of the files that were used, produced and/or modified by the run"
}

export const fileMetadata = {
    "Gaussian_Input": {
        isPlaintext: true
    },
    "Molcas_Input": {
        isPlaintext: true
    },
    "ePolyscat_Input_File": {
        isPlaintext: true
    },
    "target": {
        isPlaintext: true
    },
    isPlaintext(filename) {
        if (filename in this && "isPlaintext" in this[filename]) {
            return this[filename]["isPlaintext"];
        } else { 
            return null;
        }
    }
};

/*
export const toFromObjects = {
    fromObject(object) {
        return this[object.key].fromObject(object);
    },
    "ePolyScat_Input_File": {
        toObject(contents) {
            const originalContents = contents;

            let comments = [];

            while (/(![^\n]*$)/m.test(contents)) {
                comments.push({
                    comment: contents.match(/(![^\n]*$)/m)[0],
                    location: contents.search(/(![^\n]*$)/m)
                });

                contents = contents.replace(/(![^\n]*$)/m, "");
            }

            contents = contents
                .replace(/[ ]+$/gm, "");

            let sections = contents.replace(/\r/g, "").split(/\n[-]+[ ]*\n/gm);

            console.log(sections, contents);

            // const titleSpaces = sections[0].match(/^[ ]*|[ ]*$/, "");
            const title = sections[0].replace(/^[ ]*|[ ]*$/g, "");
            let basicData = {};
            let basicDataArray = (1 in sections)
                ? sections[1].split(/[ ]*\n[ ]*|[ ]*=[ ]/g)
                : ["LMax", "EMax", "EngForm", "VCorr"];

            for (let i = 0; i < basicDataArray.length; i += 2) {
                basicData[basicDataArray[i]] = basicDataArray[i + 1];
            }
    
            let EngForm = {};
            let EngFormDataArray = (3 in sections)
                ? sections[3].split(/[ ]*\n[ ]*|[ ]*=[ ]/g)
                : [ "iChrgMolec", "iPotFrmType" ]
            // if iPotFrmType == 3 add B.NumOrbFrm C.OrbDegn(1:NumOrbFrm) D.SymCont, SymTotal, E.OrbOccFrm(1:NumOrbFrm), F.NCoefKInt, G.CoefKInt(1:NCoefKInt)
            // if iPotFrmType == 2 add  B C OrbOccFrm(i), CoefK(i), iOrthOrb(i)
            // if iPotFrmType == 1 add B C OrbOccFrm(i), CoefK(i) 
            for (let i= 0; i < EngFormDataArray.length; i += 2) {
                  EngForm[EngFormDataArray[i]] = EngFormDataArray[i + 1];
            } 

            
            let AsyPol = {};
            let AsyPolDataArray = (5 in sections)
                ? sections[3].split(/[ ]*\n[ ]*|[ ]*=[ ]/g)
                : [ "SwitchD", "nterm" ]
            for (let i= 0; i < AsyPolDataArray.length; i += 2) {
                  AsyPol[AsyPolDataArray[i]] = AsyPolDataArray[i + 1];
            }
            // 1-nterm  A.itcen B.if itcen equal to 0 then read pcen(1:3, iterm) C.ittyp D.if ittyp equal to 1 then read apolsph(iterm) E.if ittyp equal to 2 then read apol(1, 1, iterm), apol(2, 2, iterm), apol(3, 3, iterm), apol(1, 2, iterm), apol(1, 3, iterm), apol(2, 3, iterm)

                ? sections[5].split("\n").map(ntermLine => ntermLine.split(/[ ]+/g))
                : [];

            let AsyPolDataArray2 = (5 in sections)
                ? sections[3].split(/[ ]*\n[ ]*|[ ]*=[ ]/g)
                : [ "icrtyp", "ilntyp" ]
            for (let i= 0; i < AsyPolDataArray2.length; i += 2) {
                  AsyPol[AsyPolDataArray2[i]] = AsyPolDataArray2[i + 1];
            }
            
            let ScatParams = {};
            let ScatParamsDataArray = (7 in sections)
                ? sections[3].split(/[ ]*\n[ ]*|[ ]*=[ ]/g)
                : [ "ScatEng", "FegeEng", "ScatContSym", "LMaxK" ]
            for (let i = 0; i < ScatParamsDataArray.length; i += 2) {
                ScatParams[ScatParamsDataArray[i]] = ScatParamsDataArray[i + 1];
            }

            let ProgCmds = {}; 
            let ProgCmdsDataArray = (9 in sections) 
                                ? sections[3].split(/[ ]*\n[ ]*|[ ]*=[ ]/g)
                : [ "Convert", "GetBlms", "ExpOrb", "GetPot", "Scat" ]
            for (let i = 0; i < ProgCmdsDataArray.length; i += 2) {
                ProgCmds[ProgCmdsDataArray[i]] = ProgCmdsDataArray[i + 1];
            }

            console.log(title, basicData, EngForm, AsyPol, ScatParams, ProgCmds );

            return {
                key: "ePolyscat_Input_File",
                pages: [
                    {
                      name: "Basic Data",
                      type: "table",
                      data: [
                          { 
                            name: "Title",
                            type: "text",
                            data: "title,
                            issues: (data) => (data.length == 0) ? ["Please Enter a title"] : []
                          },
                          {
                            name: "LMax",
                            type: "number",
                            data: basicData.LMax,,
                            issues: (data) => (isNaN(parseInt(data))) ? ["Enter the maximum l to be used for wave functions"] :
                                    (parseInt(data) <= 0) ? ["The number must be greater than zero"] : []
                          },
                          {
                            name: "EMax",
                            type: "number",
                            data: basicData.EMax,,
                            issues: (data) => (isNaN(parseFloat(data))) ? ["Enter the maximum asymptotic energy in eV"] : []
                          },
                            name: "EngForm",
                            type: "text",
                            data: basicData.EngForm
                            issues: (data) => (data.length == 0) ? ["Enter the charge and formula type"] : []
                          },
                          {
                            name: "VCorr,
                            type: "select",
                            data: basicData.VCorr,
                            options: ["None", "PZ", "PN", "BN", "POS-FIT"],
                          }
                      ]
                  },
                  {
                      name: "AsyPol"
                      type: ""
                      data: ""
                      columns: [


                     ]
                  },
                  {
                      name: "ScatParams"
                 
                  },
                  {
                      name: "ProgCmds"

                  }

              ],
                originalContents,
                comments,
                // targetStatesSpaces,
                // partialWavesSpaces,
                // titleSpaces,
                onchange() {

                }
            }
        },
        fromObject(object) {
            let contents =
`  ${object.pages[0].data[0].data}
   // more see for target --- but all of this commented out 
 

 
    },
    "target": {
        toObject(contents) {
            const originalContents = contents;

            let comments = [];

            while (/(![^\n]*$)/m.test(contents)) {
                comments.push({
                    comment: contents.match(/(![^\n]*$)/m)[0],
                    location: contents.search(/(![^\n]*$)/m)
                });

                contents = contents.replace(/(![^\n]*$)/m, "");
            }

            contents = contents
                .replace(/[ ]+$/gm, "");

            let sections = contents.replace(/\r/g, "").split(/\n[-]+[ ]*\n/gm);

            console.log(sections, contents);

            // const titleSpaces = sections[0].match(/^[ ]*|[ ]*$/, "");
            const title = sections[0].replace(/^[ ]*|[ ]*$/g, "");

            let basicData = {};
            let basicDataArray = (1 in sections) 
                ? sections[1].split(/[ ]*\n[ ]*|[ ]*=[ ]/g) 
                : ["coupling", "LS", "nz", ""];

            for (let i = 0; i < basicDataArray.length; i += 2) {
                basicData[basicDataArray[i]] = basicDataArray[i + 1];
            }

            let targetStates = (3 in sections) 
                ? sections[3].split("\n").map(targetStateLine => targetStateLine.split(/[ ]+/g)) 
                : [];

            // const targetStatesSpaces = sections[3]
            //     .split("\n")
            //     .map(targetStateLine => targetStateLine.match(/[ ]+/g));

            let partialWaves = (5 in sections) 
                ? sections[5].split("\n").map(partialWaveLine => partialWaveLine.split(/[ ]+/g)) 
                : [];

            // const partialWavesSpaces = sections[5]
            //     .split("\n")
            //     .map(targetStateLine => targetStateLine.match(/[ ]+/g));

            console.log(title, basicData, targetStates, partialWaves);

            return {
                key: "target",
                pages: [
                    { 
                        name: "Basic Data", 
                        type: "table", 
                        data: [
                            {
                                name: "Title",
                                type: "text",
                                data: title,
                                issues: (data) => (data.length == 0) ? ["Please Enter a title"] : []
                            },
                            { 
                                name: "Coupling Type", 
                                type: "select", 
                                data: basicData.coupling, 
                                options: ["LS", "JJ", "JK"],
                            },
                            { 
                                name: "Atomic Number",
                                type: "number",
                                data: parseInt(basicData.nz),
                                min: 1,
                                issues: (data) => 
                                    (isNaN(parseInt(data))) ? ["Enter the atomic number"] :
                                    (parseInt(data) <= 0) ? ["The atomic number must be greater than zero"] : []
                            },
                            { 
                                name: "Number of Electrons",
                                type: "number",
                                data: parseInt(basicData.nelc),
                                min: 0,
                                issues: (data) => 
                                    (isNaN(parseInt(data))) ? ["Enter the number of electrons"] :
                                    (parseInt(data) < 0) ? ["The number of electrons cannot be negative"] : []
                            }
                        ]
                    },
                    {
                        name: "Target States",
                        type: "list",
                        data: targetStates.map(([nameCFilename]) => [{ data: nameCFilename }]),
                        columns: [
                            { 
                                name: "name.c file", 
                                type: "text", 
                                issues: (data) => (data.slice(-2) != ".c") ? ["The filename must end in .c"] : []  
                            }
                        ]
                    },
                    {
                        name: "Partial Waves",
                        type: "list",
                        data: partialWaves.map(([partialWaveName, n1, n2, n3, nameCFilenameOrNo]) => [
                            { data: partialWaveName },
                            { data: parseInt(n1) },
                            { data: parseInt(n2) },
                            { data: parseInt(n3) },
                            { data: nameCFilenameOrNo }
                        ]),
                        columns: [
                            { name: "Partial Wave Name", type: "text" }, 
                            { name: "n1", type: "number" }, 
                            { name: "n2", type: "number" }, 
                            { name: "n3", type: "number" }, 
                            { name: "name.c file (or \"no\")", type: "text" }
                        ]
                    }
                ],
                originalContents,
                comments,
                // targetStatesSpaces,
                // partialWavesSpaces,
                // titleSpaces,
                onchange() {

                }
            }
        },
        fromObject(object) {
            let contents = 
`  ${object.pages[0].data[0].data}
------------------------------------------------------------------------
coupling = ${object.pages[0].data[1].data}
nz = ${object.pages[0].data[2].data}
nelc = ${object.pages[0].data[3].data}
------------------------------------------------------------------------        
ntarg = ${object.pages[1].data.length}
------------------------------------------------------------------------        
${object.pages[1].data.map(row => row.map(cell => cell.data).join("   ")).join("\n")}
------------------------------------------------------------------------
nlsp = ${object.pages[2].data.length}
------------------------------------------------------------------------        
${object.pages[2].data.map(row => row.map(cell => cell.data).join("   ")).join("\n")}
------------------------------------------------------------------------
`;
            console.log(contents);
            return contents;
        }
    },
    "Parameters": {
        toObject(parameters) {
            return {
                key: "Parameters",
                pages: [{
                    name: "Parameters",
                    type: "table",
                    data: parameters.map(parameter => ({
                        name: parameter.name,
                        type: "number",
                        data: parameter.value,
                        issues: parameter.issues,
                        dependencies: parameter.dependencies
                    }))
                }]
            }
        },
        fromObject(object) {
            return object.pages[0].data.map(item => ({
                name: item.name,
                dependencies: item.dependencies,
                value: item.data,
                issues: item.issues
            }))
        }
    }
};

*/
export const plotObjects = {
    "tr_nnn_nnn": {
        axesOptions: [
            { text: "eV", value: 0 },
            { text: "sigma", value: 1 },
            { text: "Ry", value: 2 },
            { text: "om", value: 3 },
            { text: "Ry/z^2", value: 4}
        ],
        xOptions: [
            { text: "Energy (eV)", value: 0 },
            { text: "Energy (Ry/z^2)", value: 4}
        ],
        yOptions: [
            { text: "Cross-section", value: 1 },
            { text: "Collision Strength", value: 3 },
        ]
    }
}

export const producesPlottables = ["WORKFLOW/STGF"];

const splitTarget = (contents, expectedSectionLength = 0) => {
    const result = {
        sections: contents.split(/[-]+[ ]*$/gm),
        sectionDividers: contents.match(/[-]+[ ]*$/gm)
    };

    if (result.sections.length < expectedSectionLength)
        throw new Error(`Could not parse target, expected at least ${expectedSectionLength} sections, but only found ${sections.length}`);

    return result;
}

const combineTarget = (sections, sectionDividers) => {
    return sections.map((section, i) => section + (sectionDividers[i] || "")).join("");
}

const knotDatParam = (displayName, name, isInteger) => {
    const decimalPart = (isInteger) ? "" : "(\\.\\d*)?";

    return {
        name: displayName,
        type: "number",
        min: 0,
        step: (isInteger) ? 1 : 0.00001,
        // tooltipText: (displayName != name) ? name : null,
        get: async () => {
            const contents = await store.getters["input/getContentsOfFile"]();
            let match = contents.match(new RegExp(`${name}[ ]*=[ ]*([0-9]+${decimalPart})`));
            
            if (match == null)
                match = contents.match(new RegExp(`[ ]*(\\d+${decimalPart})[ ]*==>[^(]*\\(${name}`));

            if (match != null && 1 in match)
                return match[1];
            else
                return null;
        },
        set: async (newValue) => {
            let contents = await store.getters["input/getContentsOfFile"]();
            if ((new RegExp(`${name}[ ]*=[ ]*([0-9]+${decimalPart})`)).test(contents))
                contents = contents.replace(new RegExp(`(${name}[ ]*=[ ]*[0-9]+)${decimalPart}`), `$1${newValue}`);
            else
                contents = contents.replace(new RegExp(`([ ]*)\d+${decimalPart}([ ]*==>[^(]*\(${name})`), `$1${newValue}${isInteger ? "$2" : "$3"}`);
            store.dispatch("input/setContents", { contents });
        },
        issues: (data) => 
            (isNaN(parseInt(data))) ? [`${displayName} must be a number`] :
            (parseInt(data) < 0) ? [`${displayName} cannot be negative`] : []
    };
};

export const tableObjects = {
/*
    "ePolyScat_Input_File": {
        toObject(contents) {
            const originalContents = contents;

            let comments = [];

            while (/(![^\n]*$)/m.test(contents)) {
                comments.push({
                    comment: contents.match(/(![^\n]*$)/m)[0],
                    location: contents.search(/(![^\n]*$)/m)
                });

                contents = contents.replace(/(![^\n]*$)/m, "");
            }

            contents = contents
                .replace(/[ ]+$/gm, "");

            let sections = contents.replace(/\r/g, "").split(/\n[-]+[ ]*\n/gm);

            console.log(sections, contents);

            // const titleSpaces = sections[0].match(/^[ ]*|[ ]*$/, "");
            const title = sections[0].replace(/^[ ]*|[ ]*$/g, "");
            let basicData = {};
            let basicDataArray = (1 in sections)
                ? sections[1].split(/[ ]*\n[ ]*|[ ]*=[ ]/g)
                : ["LMax", "EMax", "EngForm", "VCorr"];

            for (let i = 0; i < basicDataArray.length; i += 2) {
                basicData[basicDataArray[i]] = basicDataArray[i + 1];
            }
    
            let EngForm = {};
            let EngFormDataArray = (3 in sections)
                ? sections[3].split(/[ ]*\n[ ]*|[ ]*=[ ]/g)
                : [ "iChrgMolec", "iPotFrmType" ]
            // if iPotFrmType == 3 add B.NumOrbFrm C.OrbDegn(1:NumOrbFrm) D.SymCont, SymTotal, E.OrbOccFrm(1:NumOrbFrm), F.NCoefKInt, G.CoefKInt(1:NCoefKInt)
            // if iPotFrmType == 2 add  B C OrbOccFrm(i), CoefK(i), iOrthOrb(i)
            // if iPotFrmType == 1 add B C OrbOccFrm(i), CoefK(i) 
            for (let i= 0; i < EngFormDataArray.length; i += 2) {
                  EngForm[EngFormDataArray[i]] = EngFormDataArray[i + 1];
            } 

            
            let AsyPol = {};
            let AsyPolDataArray = (5 in sections)
                ? sections[3].split(/[ ]*\n[ ]*|[ ]*=[ ]/g)
                : [ "SwitchD", "nterm" ]
            for (let i= 0; i < AsyPolDataArray.length; i += 2) {
                  AsyPol[AsyPolDataArray[i]] = AsyPolDataArray[i + 1];
            }
            // 1-nterm  A.itcen B.if itcen equal to 0 then read pcen(1:3, iterm) C.ittyp D.if ittyp equal to 1 then read apolsph(iterm) E.if ittyp equal to 2 then read apol(1, 1, iterm), apol(2, 2, iterm), apol(3, 3, iterm), apol(1, 2, iterm), apol(1, 3, iterm), apol(2, 3, iterm)

                ? sections[5].split("\n").map(ntermLine => ntermLine.split(/[ ]+/g))
                : [];

            let AsyPolDataArray2 = (5 in sections)
                ? sections[3].split(/[ ]*\n[ ]*|[ ]*=[ ]/g)
                : [ "icrtyp", "ilntyp" ]
            for (let i= 0; i < AsyPolDataArray2.length; i += 2) {
                  AsyPol[AsyPolDataArray2[i]] = AsyPolDataArray2[i + 1];
            }
            
            let ScatParams = {};
            let ScatParamsDataArray = (7 in sections)
                ? sections[3].split(/[ ]*\n[ ]*|[ ]*=[ ]/g)
                : [ "ScatEng", "FegeEng", "ScatContSym", "LMaxK" ]
            for (let i = 0; i < ScatParamsDataArray.length; i += 2) {
                ScatParams[ScatParamsDataArray[i]] = ScatParamsDataArray[i + 1];
            }

            let ProgCmds = {}; 
            let ProgCmdsDataArray = (9 in sections) 
                                ? sections[3].split(/[ ]*\n[ ]*|[ ]*=[ ]/g)
                : [ "Convert", "GetBlms", "ExpOrb", "GetPot", "Scat" ]
            for (let i = 0; i < ProgCmdsDataArray.length; i += 2) {
                ProgCmds[ProgCmdsDataArray[i]] = ProgCmdsDataArray[i + 1];
            }

            console.log(title, basicData, EngForm, AsyPol, ScatParams, ProgCmds );
            return {
                key: 
*/
    "ePolyscat_Input_File": {
        pages: [
            {
                name: "Basic Data",
                type: "table",
                data: 
                [
                    { 
                       name: "Title",
                       type: "text",
                       get: async () => {
                          const { sections } = splitTarget(await store.getters["input/getContentsOfFile"](), 1);
                          const match = sections[0].match(/^([ ]*)([^\n\r]*[^ \n\r]|)([ ]*)$/m);

                          if (match == null)
                             return ""; 
                          else 
                             return match[2];
                       },  
                       set: async (title) => {
                          let { sections, sectionDividers } = splitTarget(await store.getters["input/getContentsOfFile"](), 1);
                          sections[0] = sections[0].replace(/^([ ]*)([^\n\r]*[^ \n\r]|)([ ]*)$/m, `$1${title}$3`);

                          const contents = combineTarget(sections, sectionDividers);

                          store.dispatch("input/setContents", { contents }); 
                       },  
                       issues: (data) => (data.length == 0) ? ["Please Enter a title"] : []
                       },
                       {
                          name: "LMax",
                          type: "number",
                          min: 1,
                          get: async () => {
                             const contents = await store.getters["input/getContentsOfFile"]();
                             const match = contents.match(/LMax[ ]+[ ]+([0-299]+)/);

                             if (match != null && 1 in match)
                                return match[1];
                             else
                                return null;
                          },
                          set: async (LMax) => {
                            let contents = await store.getters["input/getContentsOfFile"]();
                            contents = contents.replace(/LMax([ ]+)([ ]+)([0-9]+)/, `LMax$1=$2${LMax}`);
                            store.dispatch("input/setContents", { contents });
                          },
                          issues: (data) =>
                            (isNaN(parseInt(data))) ? ["Enter the maximum l to be used for wave functions"] :
                            (parseInt(data) <= 0) ? ["The number must be greater than zero"] : []
                       },
                       {
                            name: "EMax",
                            type: "number",
                            get: async () => {
                            const contents = await store.getters["input/getContentsOfFile"]();
                            const match = contents.match(/EMax[ ]+[ ]^[+-]?([0-9]*[.])?[0-9]+$/);

                            if (match != null && 1 in match)
                                return match[1];
                            else
                                return null;
                            },
                            set: async (EMax) => {
                            let contents = await store.getters["input/getContentsOfFile"]();
                            contents = contents.replace(/EMax([ ]+)([ ]+)(^[+-]?([0-9]*[.])?[0-9]+$)/, `EMax$1=$2${EMax}`);
                            store.dispatch("input/setContents", { contents });
                            },
                            issues: (data) =>
                             (isNaN(parseInt(data))) ? ["Enter the maximum asymptotic energy in eV"] : []
                       },
                       {
                            name: "EngForm",
                            type: "list",
                            async get() {
                                let contents = await store.getters["input/getContentsOfFile"]();
                                const match = contents.match(/EngForm[ ]+[ ]^[+-]?([0-9]*[.])?[0-9]+$/);

                                if (match != null && 1 in match)
                                    return match[1];
                                else
                                    return null;

                                contents = contents.replace(/\r/g, "");
                                const sections = contents.split(/[-]+[ ]*$/gm)
                                let sectionNum = 5;
            
                                if (5 in sections && /cf[ ]*=[ ]*\d+/.test(sections[5]))
                                    sectionNum = 6;

                                if (sectionNum in sections) {
                                    const rows = sections[sectionNum]
                                        .split(/[\n\r]+/g)
                                        .filter(row => !/^[ ]*$/m.test(row))
                                        .map((row, i) => ({
                                            delete: async () => {
                                                let contents = await store.getters["input/getContentsOfFile"]();
                                                const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                                                const sections = contents.split(/[-]+[ ]*$/gm);

                                                if (5 in sections) {
                                                    sections[5] = replaceNthInstance(sections[5], /[\n\r]*[^\n\r]+/m, i, "");
            
                                                    const numberRows = sections[5]
                                                        .split(/[\n\r]+/g)
                                                        .filter(row => !/^[ ]*$/m.test(row))
                                                        .length

                                                    sections[4] = sections[4].replace(/(?<=cf[ ]*=[ ]*)[^ \n\r]+/, numberRows);
                                                }
            
                                                contents = combineSeperators(sections, sectionDividers);
                                                store.dispatch("input/setContents", { contents });
                                            },
                                            cells: row
                                                .split(/[ ]+/g)
                                                .filter(row => !/^[ ]*$/m.test(row))
                                                .map((cell, j) => ({
                                                    get: async () => {
                                                        return cell;
                                                    },
                                                    set: async (data) => {
                                                        let contents = await store.getters["input/getContentsOfFile"]();

                                                        const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                                                        const sections = contents.split(/[-]+[ ]*$/gm);

                                                        if (sectionNum in sections) {
                                                            sections[sectionNum] = replaceCell(sections[sectionNum], data, i, j);
                                                        }

                                                        contents = combineSeperators(sections, sectionDividers);

                                                        store.dispatch("input/setContents", { contents });
                                                    }
                                                }))
                                        }));

                                    return rows;
                                } else {
                                    return [];
                                }
                            },
                            async addRow() {
                                let contents = await store.getters["input/getContentsOfFile"]();
                                const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                                const sections = contents.split(/[-]+[ ]*$/gm);
                                let sectionNum = 5;

                                if (5 in sections && /cf[ ]*=[ ]*\d+/.test(sections[5]))
                                    sectionNum = 6;

                                if (sectionNum in sections) {
                                    sections[sectionNum] += this.default.join("    ") + "\n";

                                    const numberRows = sections[sectionNum]
                                        .split(/[\n\r]+/g)
                                        .filter(row => !/^[ ]*$/m.test(row))
                                        .length

                                    sections[sectionNum - 1] = sections[sectionNum - 1].replace(/(?<=cf[ ]*=[ ]*)[^ \n\r]+/, numberRows);
                                }

                                contents = combineSeperators(sections, sectionDividers);
                                store.dispatch("input/setContents", { contents });
                            },
                            columns: [
                                { name: "Charge", type: "number" },
                                { name: "EngFormType", type: "number" },
                            ],
                            default: [0, 0, "no"]
                        },
                        {
                            name: "VCorr",
                            type: "select",
                            get: async () => {
                                 /*
                                 const match = contents.match(/VCorr[ ]+[ ]^[+-]?([0-9]*[.])?[0-9]+$/);

                                 if (match != null && 1 in match)
                                       return match[1];
                                 else
                                       return null;
                                 */
                            const contents = await store.getters["input/getContentsOfFile"]();
                            const match = contents.match(/VCorr[ ]+=[ ]+(None|PZ|PN|BN|POS-FIT)/);

                            if (match != null && 1 in match)
                                return match[1];
                            else
                                return null;
                            },
                            set: async (vcorrType) => {
                            let contents = await store.getters["input/getContentsOfFile"]();
                            contents = contents.replace(/vcorr([ ]+)=([ ]+)(None|PZ|PN|BN|POS-FIT)/, `vcorr$1=$2${vcorrType}`);
                            store.dispatch("input/setContents", { contents });
                            },
                            options: ["None", "PZ", "PN", "BN", "POS-FIT"],

                        }
                      ]
                    },
                    {
                      name: "AsyPol",
                      type: "list",
                      async get() {
                          let contents = await store.getters["input/getContentsOfFile"]();
                          //const match = contents.match(/[ ]AsyPol[ ]+[ ]^[+-]?([0-9]*[.])?[0-9]+$/);
                          const match = contents.match(/[ ]AsyPol[ ]/);

                            if (match != null && 1 in match)
                                return match[1];
                            else
                                return null;

                          const match2 = contents.match(/([0-9]*[.])?[0-9]+[]*# SwitchD[ ]/);

                            if (match2 != null && 1 in match2)
                                return match2[1];
                            else
                                return null;


                          contents = contents.replace(/\r/g, "");
                          const sections = contents.split(/[-]+[ ]*$/gm)
                          let sectionNum = 5;
      
                          if (5 in sections && /([0-299]+)[ ]*# nterm[ ]/.test(sections[5]))
                              sectionNum = 6;
      
                          if (sectionNum in sections) {
                              const rows = sections[sectionNum]
                                  .split(/[\n\r]+/g)
                                  .filter(row => !/^[ ]*$/m.test(row))
                                  .map((row, i) => ({
                                      delete: async () => {
                                          let contents = await store.getters["input/getContentsOfFile"]();
                                          const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                                          const sections = contents.split(/[-]+[ ]*$/gm);
      
                                          if (5 in sections) {
                                              sections[5] = replaceNthInstance(sections[5], /[\n\r]*[^\n\r]+/m, i, "");
      
                                              const numberRows = sections[5]
                                                  .split(/[\n\r]+/g)
                                                  .filter(row => !/^[ ]*$/m.test(row))
                                                  .length
      
                                              sections[4] = sections[4].replace(/(?<=lsci[ ]*=[ ]*)[^ \n\r]+/, numberRows);
                                          }
      
                                          contents = combineSeperators(sections, sectionDividers);
                                          store.dispatch("input/setContents", { contents });
                                      },
                                      cells: row
                                          .split(/[ ]+/g)
                                          .filter(row => !/^[ ]*$/m.test(row))
                                          .map((cell, j) => ({
                                              get: async () => {
                                                  return cell;
                                              },
                                              set: async (data) => {
                                                  let contents = await store.getters["input/getContentsOfFile"]();
      
                                                  const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                                                  const sections = contents.split(/[-]+[ ]*$/gm);
      
                                                  if (sectionNum in sections) {
                                                      sections[sectionNum] = replaceCell(sections[sectionNum], data, i, j);
                                                  }
      
                                                  contents = combineSeperators(sections, sectionDividers);
      
                                                  store.dispatch("input/setContents", { contents });
                                              }
                                          }))
                                  }));
      
                              return rows;
                          } else {
                              return [];
                          }
                      },
                      async addRow() {
                          let contents = await store.getters["input/getContentsOfFile"]();
                          const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                          const sections = contents.split(/[-]+[ ]*$/gm);
                          let sectionNum = 5;
      
                          if (5 in sections && /lsci[ ]*=[ ]*\d+/.test(sections[5]))
                              sectionNum = 6;
      
                          if (sectionNum in sections) {
                              sections[sectionNum] += this.default.join("    ") + "\n";
      
                              const numberRows = sections[sectionNum]
                                  .split(/[\n\r]+/g)
                                  .filter(row => !/^[ ]*$/m.test(row))
                                  .length
      
                              sections[sectionNum - 1] = sections[sectionNum - 1].replace(/(?<=lsci[ ]*=[ ]*)[^ \n\r]+/, numberRows);
                          }
      
                          contents = combineSeperators(sections, sectionDividers);
                          store.dispatch("input/setContents", { contents });
                      },
                      columns: [
                          //{ name: "Label/SwitchD", type: "text" },
                          //{ name: "nterm/Sph.Pol.", type: "number" },
                          { name: "center", type: "number" },
                          { name: "ittyp", type: "number" },
                          { name: "value", type: "number" }
                      ],
                      default: [1, 0, 0.0, "no"]
                    },                  
                    {
                      name: "ScatParams",
                      type: "list",
                      async get() {
                          const contents = await store.getters["input/getContentsOfFile"]();
                          const rowRegex = /(?<=\r?\n)[ ]*([A-Za-z0-9_]+|<[ ]*[A-Za-z0-9_]+[ ]*\|[ ]*[A-Za-z0-9_]+[ ]*>)[ ]*=[ ]*-?\d+(\.\d*)?[^#\n\r]*/g
                          const cellRegex = /<[ ]*[A-Za-z0-9_]+[ ]*\|[ ]*[A-Za-z0-9_]+[ ]*>|-?\d+(\.\d*)?|[A-Za-z0-9_]+/g;
            
                          return getGrid(contents, [rowRegex, cellRegex])
                              .map(([param, value], i) => ({
                                  async delete() {
                                      let contents = await store.getters["input/getContentsOfFile"]();
                                      // contents = getNewContentsWithoutRow(contents, removalRegex, rowRegex, i);
                                      contents = getContentsWithReplaced(contents, () => "", [[rowRegex, i]]);
                                      store.dispatch("input/setContents", { contents });
                                  },
                                  cells: [
                                      {
                                          async get() { return param; },
                                          async set(newValue) {
                                              let contents = await store.getters["input/getContentsOfFile"]();
                                              // contents = getNewContentsWithCell(contents, newValue, rowRegex, cellRegex, i, 0);
                                              contents = getContentsWithReplaced(contents, () => newValue, [[rowRegex, i], [cellRegex, 0]]);
                                              store.dispatch("input/setContents", { contents });
                                          }
                                      },
                                      {
                                          async get() { return value; },
                                          async set(newValue) {
                                              let contents = await store.getters["input/getContentsOfFile"]();
                                              // contents = getNewContentsWithCell(contents, newValue, rowRegex, cellRegex, i, 1);
                                              contents = getContentsWithReplaced(contents, () => newValue, [[rowRegex, i], [cellRegex, 1]]);
                                              store.dispatch("input/setContents", { contents });
                                          }
                                      }
                                  ]
                              }))
                      },
                      async addRow() {
                          let contents = await store.getters["input/getContentsOfFile"]();
                          contents += `\n${this.default[0]} = ${this.default[1]}`;
                                store.dispatch("input/setContents", { contents });
                      },
                      columns: [
                          {
                              name: "Parameter name",
                              type: "text",
                              issues: data =>
                                  (data.length == 0) ? ["The parameter name cannot be empty"] :
                                  (!/^[A-Za-z0-9_]+$|^<[ ]*[A-Za-z0-9_]+[ ]*\|[ ]*[A-Za-z0-9_]+[ ]*>$/.test(data)) ? ["Parameter names can only contain letters, numbers and underscores, or must be in this notation: < 2s | ks >"] : []
                          },
                          {
                              name: "Value",
                              type: "number",
                              step: "any",
                              issues: (data) => (isNaN(parseFloat(data))) ? ["Enter a number"] : []
                                  }
                              ],
                              default: ["param_name", 1.0]                       
                    },
                    {
                      name: "ProgCmds",
                      type: "list",
                      async get2() {
                          const contents = await store.getters["input/getContentsOfFile"]();
                          const sectionRegex = /(\r?\n?[^\n\r]*[^- \n\r][^\n\r]*)+\r?\n?/g; // sequence of lines that contain something that isn't a dash or space
                          const rowRegex = /\r?\n([ ]+\d+\.\d+)*[ ]*(?=\r?\n)/g
                          const cellRegex = /[^ \n\r]+/g;
                          const grid = getGrid(contents, [sectionRegex, rowRegex, cellRegex]);
      
                          if (!(4 in grid)) return [];
      
      
                          return grid
                              .map((cells, i) => ({
                                  async delete() {
                                      let contents = await store.getters["input/getContentsOfFile"]();
                                      // contents = getNewContentsWithoutRow(contents, null, rowRegex, i);
                                      contents = getContentsWithReplaced(contents, () => "", [[sectionRegex, 4], [rowRegex, i]]);
                                      store.dispatch("input/setContents", { contents });
                                  },
                                  cells: cells.map((cell, j) => ({
                                      async get() {
                                          return cell;
                                      },
                                      async set(newValue) {
                                          let contents = await store.getters["input/getContentsOfFile"]();
                                          contents = getContentsWithReplaced(contents, () => newValue, [[sectionRegex, 4], [rowRegex, i], [cellRegex, j]]);
                                          store.dispatch("input/setContents", { contents });
                                      }
                                  }))
                              }))
                      },
                      get: async () =>  {
                          let contents = await store.getters["input/getContentsOfFile"]();
      
                          contents = contents.replace(/\r/g, "");
                          const sections = contents.split(/[-]+[ ]*$/gm);
      
                          if (3 in sections) {
                              const rows = sections[3]
                                  .split(/[\n\r]+/g)
                                  .filter(row => !/^[ ]*$/m.test(row))
                                  .map((row, i) => ({
                                      delete: async () => {
                                          let contents = await store.getters["input/getContentsOfFile"]();
                                          const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                                          const sections = contents.split(/[-]+[ ]*$/gm);
      
                                          if (3 in sections) {
                                              sections[3] = replaceNthInstance(sections[3], /[\n\r]*[^\n\r]+/m, i, "");
      
                                              const numberRows = sections[3]
                                              .split(/[\n\r]+/g)
                                              .filter(row => !/^[ ]*$/m.test(row))
                                              .length
      
                                              //sections[2] = sections[2].replace(/(?<=prgcm[ ]*=[ ]*)[^ \n\r]+/, numberRows);
                                          }
      
                                          contents = combineSeperators(sections, sectionDividers);
                                          store.dispatch("input/setContents", { contents });
                                      },
                                      cells: row
                                          .split(/[ ]+/g)
                                          .filter(row => !/^[ ]*$/m.test(row))
                                          .map((cell, j) => ({
                                              get: async () => {
                                                  return cell;
                                              },
                                              set: async (data) => {
                                                  let contents = await store.getters["input/getContentsOfFile"]();
      
                                                  const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                                                  const sections = contents.split(/[-]+[ ]*$/gm);
      
                                                  if (3 in sections) {
                                                      sections[3] = replaceCell(sections[3], data, i, j);
                                                  }
      
                                                  contents = combineSeperators(sections, sectionDividers);
      
                                                  store.dispatch("input/setContents", { contents });
                                              }
                                          }))
                                  }));
      
                              return rows;
                          } else {
                              return [];
                          }
                      },
                      async addRow() {
                          let contents = await store.getters["input/getContentsOfFile"]();
                          const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                          const sections = contents.split(/[-]+[ ]*$/gm);
      
                          if (3 in sections) {
                              sections[3] += this.default.join("    ") + "\n";
      
                              const numberRows = sections[3]
                                  .split(/[\n\r]+/g)
                                  .filter(row => !/^[ ]*$/m.test(row))
                                  .length
      
                              //sections[2] = sections[2].replace(/(?<=prgcm[ ]*=[ ]*)[^ \n\r]+/, numberRows);
                          }
      
                          contents = combineSeperators(sections, sectionDividers);
                          store.dispatch("input/setContents", { contents });
                      },
                      columns: [
                          {
                              name: "Program Command",
                              type: "text",
                              issues: (data) => (data == "") ? ["This cannot be empty or delete the row"] : []
                          }
                      ],
                      default: ["GetBlms"]
                    }
        ]
    },
    "gaussian.inp": {
        pages: [

        ]
    },
    "molcas.inp": {
        pages: [

        ]
    },
    "target": {
        pages: [
            { 
                name: "Basic Data", 
                type: "table", 
                data: [
                    {
                        name: "Title",
                        type: "text",
                        get: async () => {
                            const { sections } = splitTarget(await store.getters["input/getContentsOfFile"](), 1);
                            const match = sections[0].match(/^([ ]*)([^\n\r]*[^ \n\r]|)([ ]*)$/m);

                            if (match == null)
                                return "";
                            else 
                                return match[2];
                        },
                        set: async (title) => {
                            let { sections, sectionDividers } = splitTarget(await store.getters["input/getContentsOfFile"](), 1);
                            sections[0] = sections[0].replace(/^([ ]*)([^\n\r]*[^ \n\r]|)([ ]*)$/m, `$1${title}$3`);

                            const contents = combineTarget(sections, sectionDividers);

                            store.dispatch("input/setContents", { contents });
                        },
                        issues: (data) => (data.length == 0) ? ["Please Enter a title"] : []
                    },
                    {
                        name: "Coupling Type", 
                        type: "select", 
                        get: async () => {
                            const contents = await store.getters["input/getContentsOfFile"]();
                            const match = contents.match(/coupling[ ]+=[ ]+(LS|JJ|JK)/);
                            
                            if (match != null && 1 in match)
                                return match[1];
                            else
                                return null;
                        },
                        set: async (couplingType) => {
                            let contents = await store.getters["input/getContentsOfFile"]();
                            contents = contents.replace(/coupling([ ]+)=([ ]+)(LS|JJ|JK)/, `coupling$1=$2${couplingType}`);
                            store.dispatch("input/setContents", { contents });
                        },
                        options: ["LS", "JJ", "JK"],
                    },
                    {
                        name: "Atomic Number",
                        type: "number",
                        min: 1,
                        get: async () => {
                            const contents = await store.getters["input/getContentsOfFile"]();
                            const match = contents.match(/nz[ ]+=[ ]+([0-9]+)/);
                            
                            if (match != null && 1 in match)
                                return match[1];
                            else
                                return null;
                        },
                        set: async (nz) => {
                            let contents = await store.getters["input/getContentsOfFile"]();
                            contents = contents.replace(/nz([ ]+)=([ ]+)([0-9]+)/, `nz$1=$2${nz}`);
                            store.dispatch("input/setContents", { contents });
                        },
                        issues: (data) => 
                            (isNaN(parseInt(data))) ? ["Enter the atomic number"] :
                            (parseInt(data) <= 0) ? ["The atomic number must be greater than zero"] : []
                    },
                    { 
                        name: "Number of Electrons",
                        type: "number",
                        min: 0,
                        get: async () => {
                            const contents = await store.getters["input/getContentsOfFile"]();
                            const match = contents.match(/nelc[ ]+=[ ]+([0-9]+)/);
                            
                            if (match != null && 1 in match)
                                return match[1];
                            else
                                return null;
                        },
                        set: async (nelc) => {
                            let contents = await store.getters["input/getContentsOfFile"]();
                            contents = contents.replace(/nelc([ ]+)=([ ]+)([0-9]+)/, `nelc$1=$2${nelc}`);
                            store.dispatch("input/setContents", { contents });
                        },
                        issues: (data) => 
                            (isNaN(parseInt(data))) ? ["Enter the number of electrons"] :
                            (parseInt(data) < 0) ? ["The number of electrons cannot be negative"] : []
                    }
                ]
            },
            {
                name: "Target States",
                type: "list",
                async get2() {
                    const contents = await store.getters["input/getContentsOfFile"]();
                    const sectionRegex = /(\r?\n?[^\n\r]*[^- \n\r][^\n\r]*)+\r?\n?/g; // sequence of lines that contain something that isn't a dash or space
                    const rowRegex = /\r?\n([ ]+\d+\.\d+)*[ ]*(?=\r?\n)/g
                    const cellRegex = /[^ \n\r]+/g;
                    const grid = getGrid(contents, [sectionRegex, rowRegex, cellRegex]);

                    if (!(4 in grid)) return [];

                    
                    return grid
                        .map((cells, i) => ({
                            async delete() {
                                let contents = await store.getters["input/getContentsOfFile"]();
                                // contents = getNewContentsWithoutRow(contents, null, rowRegex, i);
                                contents = getContentsWithReplaced(contents, () => "", [[sectionRegex, 4], [rowRegex, i]]);
                                store.dispatch("input/setContents", { contents });
                            },
                            cells: cells.map((cell, j) => ({
                                async get() { 
                                    return cell;
                                },
                                async set(newValue) {
                                    let contents = await store.getters["input/getContentsOfFile"]();
                                    contents = getContentsWithReplaced(contents, () => newValue, [[sectionRegex, 4], [rowRegex, i], [cellRegex, j]]);
                                    store.dispatch("input/setContents", { contents });
                                }
                            }))
                        }))
                },
                get: async () =>  {
                    let contents = await store.getters["input/getContentsOfFile"]();

                    contents = contents.replace(/\r/g, "");
                    const sections = contents.split(/[-]+[ ]*$/gm);
                    
                    if (3 in sections) {
                        const rows = sections[3]
                            .split(/[\n\r]+/g)
                            .filter(row => !/^[ ]*$/m.test(row))
                            .map((row, i) => ({
                                delete: async () => {
                                    let contents = await store.getters["input/getContentsOfFile"]();
                                    const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                                    const sections = contents.split(/[-]+[ ]*$/gm);

                                    if (3 in sections) {
                                        sections[3] = replaceNthInstance(sections[3], /[\n\r]*[^\n\r]+/m, i, "");

                                        const numberRows = sections[3]
                                        .split(/[\n\r]+/g)
                                        .filter(row => !/^[ ]*$/m.test(row))
                                        .length

                                        sections[2] = sections[2].replace(/(?<=ntarg[ ]*=[ ]*)[^ \n\r]+/, numberRows);
                                    }

                                    contents = combineSeperators(sections, sectionDividers);
                                    store.dispatch("input/setContents", { contents });
                                },
                                cells: row
                                    .split(/[ ]+/g)
                                    .filter(row => !/^[ ]*$/m.test(row))
                                    .map((cell, j) => ({
                                        get: async () => {
                                            return cell;
                                        },
                                        set: async (data) => {
                                            let contents = await store.getters["input/getContentsOfFile"]();
                        
                                            const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                                            const sections = contents.split(/[-]+[ ]*$/gm);
                        
                                            if (3 in sections) {
                                                sections[3] = replaceCell(sections[3], data, i, j);
                                            }
                        
                                            contents = combineSeperators(sections, sectionDividers);
                        
                                            store.dispatch("input/setContents", { contents });
                                        }
                                    }))
                            }));

                        return rows;
                    } else {
                        return [];
                    }
                },
                async addRow() {
                    let contents = await store.getters["input/getContentsOfFile"]();
                    const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                    const sections = contents.split(/[-]+[ ]*$/gm);
                    
                    if (3 in sections) {
                        sections[3] += this.default.join("    ") + "\n";

                        const numberRows = sections[3]
                            .split(/[\n\r]+/g)
                            .filter(row => !/^[ ]*$/m.test(row))
                            .length

                        sections[2] = sections[2].replace(/(?<=ntarg[ ]*=[ ]*)[^ \n\r]+/, numberRows);
                    }

                    contents = combineSeperators(sections, sectionDividers);
                    store.dispatch("input/setContents", { contents });
                },
                columns: [
                    { 
                        name: "name.c file", 
                        type: "text",
                        issues: (data) => (data == "") ? ["The filename cannot be empty"] : []  
                    }
                ],
                default: ["filename.c"]
            },
            {
                name: "Partial Waves",
                type: "list",
                async get() {
                    let contents = await store.getters["input/getContentsOfFile"]();
                    contents = contents.replace(/\r/g, "");
                    const sections = contents.split(/[-]+[ ]*$/gm)
                    let sectionNum = 5;

                    if (5 in sections && /nlsp[ ]*=[ ]*\d+/.test(sections[5]))
                        sectionNum = 6;

                    if (sectionNum in sections) {
                        const rows = sections[sectionNum]
                            .split(/[\n\r]+/g)
                            .filter(row => !/^[ ]*$/m.test(row))
                            .map((row, i) => ({
                                delete: async () => {
                                    let contents = await store.getters["input/getContentsOfFile"]();
                                    const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                                    const sections = contents.split(/[-]+[ ]*$/gm);

                                    if (5 in sections) {
                                        sections[5] = replaceNthInstance(sections[5], /[\n\r]*[^\n\r]+/m, i, "");

                                        const numberRows = sections[5]
                                            .split(/[\n\r]+/g)
                                            .filter(row => !/^[ ]*$/m.test(row))
                                            .length
                
                                        sections[4] = sections[4].replace(/(?<=nlsp[ ]*=[ ]*)[^ \n\r]+/, numberRows);
                                    }

                                    contents = combineSeperators(sections, sectionDividers);
                                    store.dispatch("input/setContents", { contents });
                                },
                                cells: row
                                    .split(/[ ]+/g)
                                    .filter(row => !/^[ ]*$/m.test(row))
                                    .map((cell, j) => ({
                                        get: async () => {
                                            return cell;
                                        },
                                        set: async (data) => {
                                            let contents = await store.getters["input/getContentsOfFile"]();
                        
                                            const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                                            const sections = contents.split(/[-]+[ ]*$/gm);
                        
                                            if (sectionNum in sections) {
                                                sections[sectionNum] = replaceCell(sections[sectionNum], data, i, j);
                                            }
                        
                                            contents = combineSeperators(sections, sectionDividers);
                        
                                            store.dispatch("input/setContents", { contents });
                                        }
                                    }))
                            }));

                        return rows;
                    } else {
                        return [];
                    }
                },
                async addRow() {
                    let contents = await store.getters["input/getContentsOfFile"]();
                    const sectionDividers = contents.match(/[-]+[ ]*$/gm);
                    const sections = contents.split(/[-]+[ ]*$/gm);
                    let sectionNum = 5;

                    if (5 in sections && /nlsp[ ]*=[ ]*\d+/.test(sections[5]))
                        sectionNum = 6;

                    if (sectionNum in sections) {
                        sections[sectionNum] += this.default.join("    ") + "\n";

                        const numberRows = sections[sectionNum]
                            .split(/[\n\r]+/g)
                            .filter(row => !/^[ ]*$/m.test(row))
                            .length

                        sections[sectionNum - 1] = sections[sectionNum - 1].replace(/(?<=nlsp[ ]*=[ ]*)[^ \n\r]+/, numberRows);
                    }

                    contents = combineSeperators(sections, sectionDividers);
                    store.dispatch("input/setContents", { contents });
                },
                columns: [
                    { name: "Label", type: "text" }, 
                    { name: "L", type: "number" }, 
                    { name: "S", type: "number" }, 
                    { name: "parity", type: "number" }
                ],
                default: ["label", 0, 0, 0, "no"]
            }
        ]
    },
    "knot.dat": {
        pages: [
            {
                name: "Options",
                type: "table",
                data: [
                    knotDatParam("Grid type", "grid_type", true),
                    knotDatParam("Order of splines", "ks", true),
                    knotDatParam("Number of splines", "ns", true),
                    knotDatParam("Nuclear charge", "z", false),
                    knotDatParam("Step size from 0 to 1, for z*r", "h", false),
                    knotDatParam("Maximum step size for r", "hmax", false),
                    knotDatParam("Maximum r", "rmax", false),
                ]
            },
            {
                name: "Knot data",
                type: "list",
                contiguousCells: true,
                async get() {
                    const contents = await store.getters["input/getContentsOfFile"]();
                    const gridRegex = /(?<=\*\*\*)[^]+(?=\*\*\*)/g;
                    const rowRegex = /\r?\n([ ]+\d+\.\d+)*[ ]*(?=\r?\n)/g;
                    const cellRegex = /[ ]*\d+(\.\d*)?/g;
                    const grid = getGrid(contents, [gridRegex, rowRegex, cellRegex])[0];

                    return grid
                        .map((cells, i) => ({
                            async delete() {
                                let contents = await store.getters["input/getContentsOfFile"]();
                                // contents = getNewContentsWithoutRow(contents, null, rowRegex, i);
                                contents = getContentsWithReplaced(contents, () => "", [[gridRegex, 0], [rowRegex, i]]);
                                store.dispatch("input/setContents", { contents });
                            },
                            cells: cells.map((cell, j) => ({
                                async get() { 
                                    return parseFloat(cell.replace(" ", ""));
                                },
                                async set(newValue) {
                                    newValue = parseFloat(newValue);
                                    const numSpaces = 5 - ((newValue < 1) ? 0 : Math.floor(Math.log10(Math.floor(newValue))))
                                    newValue = " ".repeat(numSpaces) + newValue.toFixed(5);
                                    let contents = await store.getters["input/getContentsOfFile"]();
                                    contents = getContentsWithReplaced(contents, () => newValue, [[gridRegex, 0], [rowRegex, i], [cellRegex, j]]);
                                    store.dispatch("input/setContents", { contents });
                                }
                            })),
                            addCell: async () => {
                                let contents = await store.getters["input/getContentsOfFile"]();
                                contents = getContentsWithReplaced(contents, (row) => `${row}     0.00000`, [[gridRegex, 0], [rowRegex, i]]);

                                if (cells.length >= 4)
                                    contents = getContentsWithReplaced(contents, grid => `${grid}\n `, [[gridRegex, 0]]);

                                store.dispatch("input/setContents", { contents });
                            },
                            async removeCell() {
                                let contents = await store.getters["input/getContentsOfFile"]();

                                if (cells.length == 0)
                                    contents = getContentsWithReplaced(contents, () => "", [[gridRegex, 0], [rowRegex, i]]);

                                if (cells.length > 0) {
                                    contents = getContentsWithReplaced(contents, (_cell) => ``, [
                                        [gridRegex, 0], [rowRegex, i], [cellRegex, cells.length - 1]
                                    ]);
                                } else if (grid.length > 0) {
                                    contents = getContentsWithReplaced(contents, (_cell) => ``, [
                                        [gridRegex, 0], [rowRegex, i - 1], [cellRegex, grid[i - 1].length - 1]
                                    ]);
                                }

                                store.dispatch("input/setContents", { contents });
                            }
                        }))
                },
                // async addRow() {
                //     let contents = await store.getters["input/getContentsOfFile"]();
                //     const gridRegex = /(?<=\*\*\*)[^]+(?=\*\*\*)/g;
                //     const rowRegex = /\r?\n([ ]+\d+\.\d+)*[ ]*(?=\r?\n)/g;
                //     const cellRegex = /[ ]*\d+(\.\d*)?/g;
                    
                //     for (let i = 0; i < 5; i++) {
                //         let grid = getGrid(contents, [gridRegex, rowRegex, cellRegex])[0];

                //         contents = getContentsWithReplaced(contents, (row) => `${row}     0.00000`, [[gridRegex, 0], [rowRegex, grid.length - 1]]);

                //         if (grid[grid.length - 1].length >= 4)
                //             contents = getContentsWithReplaced(contents, grid => `${grid}\n `, [[gridRegex, 0]]);
                //     }

                //     store.dispatch("input/setContents", { contents });
                //     // let contents = await store.getters["input/getContentsOfFile"]();
                //     // const gridRegex = /(?<=\*\*\*)[^]+(?=\*\*\*)/g;
                //     // // contents += `\n${this.default.join("")}`;
                //     // contents = getContentsWithReplaced(contents, grid => `${grid}\n `, [[gridRegex, 0]]);
                //     // store.dispatch("input/setContents", { contents });
                // },
                columns: [
                    { 
                        name: "", 
                        type: "number", 
                        step: 0.00001,
                        issues: (data) => (isNaN(parseFloat(data))) ? ["Enter a number"] : []
                    },
                    { 
                        name: "", 
                        type: "number", 
                        step: 0.00001,
                        issues: (data) => (isNaN(parseFloat(data))) ? ["Enter a number"] : []
                    },
                    { 
                        name: "", 
                        type: "number", 
                        step: 0.00001,
                        issues: (data) => (isNaN(parseFloat(data))) ? ["Enter a number"] : []
                    },
                    { 
                        name: "", 
                        type: "number", 
                        step: 0.00001,
                        issues: (data) => (isNaN(parseFloat(data))) ? ["Enter a number"] : []
                    },
                    { 
                        name: "", 
                        type: "number", 
                        step: 0.00001,
                        issues: (data) => (isNaN(parseFloat(data))) ? ["Enter a number"] : []
                    },
                ],
                default: ["     0.00000", "     0.00000", "     0.00000", "     0.00000", "     0.00000"]
            }
        ]
    },
    "bsr_par": {
        pages: [
            {
                name: "Extra parameters",
                type: "list",
                async get() {
                    const contents = await store.getters["input/getContentsOfFile"]();
                    const rowRegex = /(?<=\r?\n)[ ]*([A-Za-z0-9_]+|<[ ]*[A-Za-z0-9_]+[ ]*\|[ ]*[A-Za-z0-9_]+[ ]*>)[ ]*=[ ]*-?\d+(\.\d*)?[^#\n\r]*/g
                    const cellRegex = /<[ ]*[A-Za-z0-9_]+[ ]*\|[ ]*[A-Za-z0-9_]+[ ]*>|-?\d+(\.\d*)?|[A-Za-z0-9_]+/g;

                    return getGrid(contents, [rowRegex, cellRegex])
                        .map(([param, value], i) => ({
                            async delete() {
                                let contents = await store.getters["input/getContentsOfFile"]();
                                // contents = getNewContentsWithoutRow(contents, removalRegex, rowRegex, i);
                                contents = getContentsWithReplaced(contents, () => "", [[rowRegex, i]]);
                                store.dispatch("input/setContents", { contents });
                            },
                            cells: [
                                {
                                    async get() { return param; },
                                    async set(newValue) {
                                        let contents = await store.getters["input/getContentsOfFile"]();
                                        // contents = getNewContentsWithCell(contents, newValue, rowRegex, cellRegex, i, 0);
                                        contents = getContentsWithReplaced(contents, () => newValue, [[rowRegex, i], [cellRegex, 0]]);
                                        store.dispatch("input/setContents", { contents });
                                    }
                                },
                                {
                                    async get() { return value; },
                                    async set(newValue) {
                                        let contents = await store.getters["input/getContentsOfFile"]();
                                        // contents = getNewContentsWithCell(contents, newValue, rowRegex, cellRegex, i, 1);
                                        contents = getContentsWithReplaced(contents, () => newValue, [[rowRegex, i], [cellRegex, 1]]);
                                        store.dispatch("input/setContents", { contents });
                                    }
                                }
                            ]
                        }))
                },
                async addRow() {
                    let contents = await store.getters["input/getContentsOfFile"]();
                    contents += `\n${this.default[0]} = ${this.default[1]}`;
                    store.dispatch("input/setContents", { contents });
                },
                columns: [
                    { 
                        name: "Parameter name", 
                        type: "text", 
                        issues: data => 
                            (data.length == 0) ? ["The parameter name cannot be empty"] : 
                            (!/^[A-Za-z0-9_]+$|^<[ ]*[A-Za-z0-9_]+[ ]*\|[ ]*[A-Za-z0-9_]+[ ]*>$/.test(data)) ? ["Parameter names can only contain letters, numbers and underscores, or must be in this notation: < 2s | ks >"] : []
                    },
                    { 
                        name: "Value", 
                        type: "number", 
                        step: "any",
                        issues: (data) => (isNaN(parseFloat(data))) ? ["Enter a number"] : []
                    }
                ],
                default: ["param_name", 1.0]
            }
        ]
    },
    "Parameters": {
        key: "Parameters",
        pages: [{
            name: "Parameters",
            type: "table",
            get: async () => {
                return store
                    .getters["input/getInputs"]().parameters
                    .map(parameter => {
                        return {
                            name: parameter.description || parameter.name,
                            type: "number",
                            get: () => parameter.value,
                            set: (value) => store.commit("input/SET_PARAMETER", {
                                name: parameter.name,
                                value: value
                            }),
                            issues: parameter.issues,
                            ...parametersData(parameter)
                        }
                    })
            }
        }]
    }
}

const parametersData = (parameter) => ({
    "cross section units": {
        name: "Cross section units", 
        type: "select", 
        get: () => parameter.value || 0,
        set: (value) => store.commit("input/SET_PARAMETER", {
            name: parameter.name,
            value: value
        }),
        options: [
            {text: "Atomic Units", value: 0},
            {text: "10^-16 cm^2", value: 1},
            {text: "10^-18 cm^2", value: 2}
        ],
        tooltip: "Select the units to be used for the cross sections"
    },
    "ii1": { min: 1, tooltip: "Enter the first initial scattering state" },
    "ii2": { min: 1, tooltip: "Enter the second initial scattering state" },
    "ff1": { min: 1, tooltip: "Enter the first final scattering state" },
    "ff2": { min: 1, tooltip: "Enter the second final scattering state" },
    "klsp1=": { min: 1, tooltip: "Enter the first index of the partial wave to have its state calculated" },
    "klsp2=": { min: 1, tooltip: "Enter the last index of the partial wave to have its state calculated" }
}[parameter.name]);

const getGrid = (contents, regexes) => {
    if (regexes.length == 0)
        return contents
    else 
        return [...contents.matchAll(regexes[0])].map(([match]) => {
            return getGrid(match, regexes.slice(1));
        })
}

const getContentsWithReplaced = (contents, replacer, regexes) => {
    if (regexes.length == 0) {
        return replacer(contents);
    } else {
        const [regex, index] = regexes[0];
        const match = [...contents.matchAll(regex)][index];
        const contentsA = contents.slice(0, match.index);
        const contentsB = contents.slice(match.index + match[0].length);

        return contentsA + getContentsWithReplaced(match[0], replacer, regexes.slice(1)) + contentsB;
    }
}

// const getGrid = (
//     contents, 
//     rowRegex, cellRegex
// ) => {
//     const cells = [...contents.matchAll(rowRegex)].map(([row]) => {
//         return [...row.matchAll(cellRegex)].map(([cell]) => cell);
//     });

//     return cells;
// }

const getNewContentsWithoutRow = (
    contents, deleteSpecificRegex = null,
    rowRegex,
    rowNumber
) => {
    const row = [...contents.matchAll(rowRegex)][rowNumber];
    const contentsA = contents.slice(0, row.index);
    const contentsB = contents.slice(row.index + row[0].length);
    let remainingLine = "";

    if (deleteSpecificRegex)
        remainingLine = row[0].replace(deleteSpecificRegex, "");

    return contentsA + remainingLine + contentsB
}

const getNewContentsWithCell = (
    contents, newValue,
    rowRegex, cellRegex, 
    rowNumber, columnNumber
) => {
    const row = [...contents.matchAll(rowRegex)][rowNumber];
    const contentsA = contents.slice(0, row.index)
    const contentsB = contents.slice(row.index + row[0].length);
    const cell = [...row[0].matchAll(cellRegex)][columnNumber];
    const rowA = row[0].slice(0, cell.index);
    const rowB = row[0].slice(cell.index + cell[0].length);

    return contentsA + rowA + newValue + rowB + contentsB;
}

const getNewContentsWithReplacedRow = (
    contents, rowReplacer = row => row,
    rowRegex,
    rowNumber
) => {
    const row = [...contents.matchAll(rowRegex)][rowNumber];
    const contentsA = contents.slice(0, row.index)
    const contentsB = contents.slice(row.index + row[0].length);

    return contentsA + rowReplacer(row) + contentsB;
}

const replaceCell = (string, replaceString, rowNumber, columnNumber) => {
    const lineSeperators = string.match(/[\n\r]+/g) || [];
    const lines = string.split(/[\n\r]+/g);
    let lineNumber = 0;

    for (let i = 0; lineNumber <= rowNumber; i++)
        if (!/^[ ]*$/m.test(lines[i]))
            lineNumber++;

    const cellSeperators = lines[lineNumber].match(/[ ]+/g) || [];
    const cells = lines[lineNumber].split(/[ ]+/g);

    cells[columnNumber] = replaceString;

    lines[lineNumber] = combineSeperators(cells, cellSeperators);

    return combineSeperators(lines, lineSeperators);
}

const combineSeperators = (seperated, seperators) => {
    seperators.push("");

    return seperated.reduce((acc, str, i) => acc + str + seperators[i], "");
}

const replaceNthInstance = (string, regex, n, replaceString) => {
    let result = "";
    let i = 0;

    while (i < n) {
        let index = string.search(regex);

        result += string.slice(0, index) + string.match(regex)[0];
        string = string.slice(index).replace(regex, "");

        i++;
    }

    return result + string.replace(regex, replaceString);
}
