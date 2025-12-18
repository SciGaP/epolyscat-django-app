import { InputService } from "@/service/epolyscat-service";
import { descriptions, tableObjects, plotObjects } from "@/fileData";

const state = {
    pathLabels: [],
    inputFiles: {},
    outputFiles: {},
    parameters: [],
    path: []
}

const actions = {
    // path should already be set to its default value by the time this is executed
    async fetchPathLabels({commit, state}) {
        let applicationInputs = await InputService.fetchApplicationInputs();
        const parameter_dependencies = applicationInputs
            .filter(input => !("config" in input.metaData.editor) && input.type <= 1)
            .reduce((dependencies, parameter) => dependencies.concat(parameter.metaData.editor.dependencies.show.OR), []);

        let inputFiles = {};
        let pathLabels = [{}, {}, {}, {}];
        let labelsKeyIndeces = {};
        let parameters = [];

        // applicationInputs = applicationInputs.filter(input => "config" in input.metaData.editor || input.type > 1);

        // Place the parameters input in between the top level options and the actual parameters
        applicationInputs.splice(10, 0, {name: "Parameters", type: 1, metaData: { editor: { dependencies: { show: { OR: parameter_dependencies } } } } })

        for (let input of applicationInputs) {
            console.log('input-storage.store:31: Current input is', input)
            let dependencies = [];
            let pathLabelsIndex = 0;

            // if (input.type < 2 && !("config" in input.metaData.editor) && input.name != "Parameters") 
            //     input.metaData.editor.dependencies.show.OR = [{"Parameters": { value: "Parameters" }}];
            
            if ("dependencies" in input.metaData.editor) {
                for (const dependencyObj of input.metaData.editor.dependencies.show.OR) {
                    const dependency = Object.values(dependencyObj)[0].value;
                    const dependencyKey = Object.keys(dependencyObj)[0];
                    pathLabelsIndex = labelsKeyIndeces[dependencyKey] + 1;

                    pathLabels[pathLabelsIndex - 1][dependency].forEach(dependencyPath => dependencies.push([...dependencyPath, dependency]));
                }
            }

            labelsKeyIndeces[input.name] = pathLabelsIndex;

            if (dependencies.length == 0)
                dependencies = [[]];
            
            if ("config" in input.metaData.editor) {
                input.metaData.editor.config.options.forEach(({ value }) => 
                    pathLabels[pathLabelsIndex][value] = dependencies
                );
            } else if (input.type > 1 || input.name == "Parameters") {
                pathLabels[pathLabelsIndex][input.name] = dependencies;
            } else {
                parameters.push({
                    name: input.name,
                    dependencies,
                    value: input.value,
                    description: input.userFriendlyDescription,
                    issues: (data) => 
                        (isNaN(parseInt(data))) ? [`Enter a value for ${input.name}`] :
                        (parseInt(data) != data) ? [`${input.name} must be an integer`] : []
                });
            }

            if (input.type == 3 || input.type == 4) {
                if ("userFriendlyDescription" in input && input["userFriendlyDescription"] != null)
                    descriptions[input.name] = input.userFriendlyDescription;

                inputFiles[input.name] = {
                    isMultiFileInput: input.type == 4,
                    files: [],
                    isValid: function(_inputFiles) { 
                        return this.files.filter(file => !file.deleted).length > 0 
                        // && this.files.every(file => file.)
                    },
                    dependencies
                };
            }
        }

        // pathLabels = pathLabels.slice(0, 2);
        commit("SET_INPUTS", { pathLabels, inputFiles, parameters, outputFiles: {} });
    },
    addFiles({commit}, { files }) {
        for (const file of files) {
            const inputFileName = getProperFilenameOf(file.name, state.inputFiles, state.path);

            if (inputFileName != null)
                commit("ADD_TO_INPUT_FILE", { file, inputFileName })
        }
    },
    setContents({commit, getters}, { contents }) {
        const file = getters["getFile"];

        file.contents = contents;
        file.unchanged = false;
    },
    async loadOutputsIntoInputs({commit, getters}, { runId }) {
        const outputFileObjs = await InputService.fetchOutputs(runId);
        const inputFiles = getters["getInputs"]({customPath: []}).inputFiles;

        for (const outputFileObj of outputFileObjs) {
            if (outputFileObj != null && typeof outputFileObj == "object" && !("dataProductURI" in outputFileObj))
                outputFileObj.dataProductURI = outputFileObj["data-product-uri"]

            const inputFileName = getProperFilenameOf(outputFileObj.name, inputFiles, []);

            if (inputFileName) commit("ADD_TO_INPUT_FILE", {
                file: outputFileObj,
                inputFileName
            })
        }

        // run.inputs[""]
    },
}

const mutations = {
    SET_INPUTS(state, { 
        inputFiles = state.inputFiles, 
        pathLabels = state.pathLabels, 
        parameters = state.parameters,
        outputFiles = state.outputFiles
    } = state) {
        state.inputFiles = inputFiles;
        state.pathLabels = pathLabels;
        state.parameters = parameters;
        state.outputFiles = outputFiles;
    },
    ADD_TO_INPUT_FILE(state, { file, inputFileName }) {
        // console.log("adding input file:", file);
        const inputFile = state.inputFiles[inputFileName];

        const indexToBeReplaced = inputFile.files
            .findIndex(other_file => other_file.name == file.name);

        if (indexToBeReplaced == -1 && !inputFile.isMultiFileInput && file.replaceCurrent) {
            // Deal with cases where the new file you're uploading has a different name, but only one can exist
            inputFile.files[0].deleted = true;
            inputFile.files.unshift(file);
        } else if (indexToBeReplaced == -1) {
            inputFile.files.push(file);
            inputFile.files.sort((a,b) => a.name.localeCompare(b.name));
        } else if (file.replaceCurrent || inputFile.files[indexToBeReplaced].deleted)
            inputFile.files.splice(indexToBeReplaced, 1, file);
    },
    // Assumes filename and properFilename are in pathLabels and inputFiles
    REMOVE_FILE(state, { filename, inputFileName }) {
        // delete state.pathLabels[3][filename];

        const file = state.inputFiles[inputFileName].files
            .find(file => file.name == filename);

        file.deleted = true;
        file.unchanged = false;
    },
    SET_PARAMETER(state, { name, value }) {
        state.parameters
            .find(parameter => parameter.name == name).value = value;
    },
    SET_PATH(state, { path }) {
        // console.log("setting path to", path);
        state.path = path;
    },
}

const getters = {
    getInputs: (state) => {
        // filters the inputs by the path
        return ({removeDeleted = true, customPath = null} = {}) => {
            const path = (customPath != null) ? customPath : state.path;

            const inputFiles = Object
                .entries(state.inputFiles)
                .filter(([_, inputFile]) => dependenciesMatchesPath(inputFile.dependencies, path))
                .reduce((filteredInputFiles, [filename, inputFile]) => {
                    filteredInputFiles[filename] = {
                        ...inputFile,
                        files: inputFile.files.filter(file => !file.deleted || !removeDeleted)
                    };
                    return filteredInputFiles;
                }, {})

            const parameters = state.parameters
                .filter(parameter => dependenciesMatchesPath(parameter.dependencies, path));

            const pathLabels = state.pathLabels.map(pathLabel => {
                return Object
                    .entries(pathLabel)
                    .filter(([_, dependencies]) => dependenciesMatchesPath(dependencies, path))
                    .map(([pathOption, _]) => pathOption);
            });
            
            return { inputFiles, pathLabels, parameters, outputFiles: state.outputFiles };
        }
    },
    getPreparedInputs: (state, getters) => {
        return async ({ prepareForCreation } = {}) => {
            const {inputFiles, parameters} = getters["getInputs"]({ 
                removeDeleted: false 
            });
    
            const preparedInputFiles = await Promise.all(Object.entries(inputFiles)
                .map(async ([inputFileName, inputFile]) => {
                    const files = await Promise.all(inputFile.files
                        .filter(file => !prepareForCreation || !file.deleted)
                        .filter(file => prepareForCreation || !file.unchanged)
                        .map(async file => {
                            const sendContents = file.isFromComputer || !!file.contents;
                            let contents = null;
    
                            if (sendContents)
                                contents = await InputService.fetchFileContents(file);
    
                            const isPlaintext = contents != null;
    
                            if (sendContents && !isPlaintext)
                                contents = await InputService.fetchFileContentsBase64(file);
    
                            return {
                                name: file.name,
                                deleted: file.deleted || false,
                                dataProductURI: file.dataProductURI,
                                contents,
                                isPlaintext
                            };
                        }
                    ))
    
                    return {
                        type: "files",
                        name: inputFileName,
                        files
                    }
                })
            );

            // console.log(preparedInputFiles);
            
            const preparedParameters = parameters.map(parameter => {
                return { 
                    type: "parameter", name: parameter.name, 
                    value: parameter.value, issues: parameter.issues 
                };
            })
            
            const preparedRunType = [
                { type: "radio input", name: "Calculation_Type", value: state.path[0] },
                (state.path[0] == "MODULE")
                    ? { type: "radio input", name: "EPOLYSCAT_Application_Module", value: state.path[1] }
                : (state.path[0] == "UTILITY")
                    ? { type: "radio input", name: "Application_Utility", value: state.path[1] }
                    : { type: "radio input", name: "Application_Workflow", value: state.path[1] }
            ];

            // console.log(preparedInputFiles.concat(preparedParameters, preparedRunType));
    
            return preparedInputFiles.concat(preparedParameters, preparedRunType);
        }
    },
    getFile: (state) => {
        const files = {...state.inputFiles, ...state.outputFiles};
        let file = null;

        if (state.path.length == 4 && Array.isArray(files[state.path[2]].files))
            file = files[state.path[2]].files
                .find(file => file.name == state.path[3]);
        else if (state.path.length == 5)
            file = files[state.path[2]]
                .files[state.path[3]]
                .files.find(file => file.name == state.path[4]);
        // console.log("File gotten: ", file)
        return file;
    },
    getFilesInCategory: (state) => {
        return (properFilename) => {
            const files = {...state.inputFiles, ...state.outputFiles};

            if (properFilename in files) {
                return files[properFilename].files;
            } else {
                for (const category in files)
                    if ("files" in files[category] && properFilename in files[category].files)
                        return files[category].files[properFilename].files;
            }
        }
    },
    getContentsOfFile: (state, getters) => {
        return async () => {
            const file = getters["getFile"];
            // console.log(await InputService.fetchFileContents(file));
            // console.log("getting contents of", file)
    
            return await InputService.fetchFileContents(file);
        }
    },
    checkIfDuplicate: (state) => {
        return (filename) => {
            const properFilename = getProperFilenameOf(filename, state.inputFiles, state.path);

            return properFilename != null && 
                state.inputFiles[properFilename].files.some(file => 
                    !file.deleted && (
                        file.name == filename
                    )
                );
        };
    },
    multipleFilesCanBeUploaded: (state) => {
        return () => {
            if (state.path[state.path.length - 1] == "Upload") {

                return state.inputFiles[state.path[state.path.length - 2]].isMultiFileInput;
            } else {
                return false;
            }
        }
    },
    getPath: (state) => state.path,
    getSelectedInfo: (state, getters) => {
        return () => {
            const file = getters["getFile"];
            const selected = state.path[state.path.length - 1];
            const properFilename = getters["getProperFilenameOf"](selected)
            let tableObject = null;
            let plotObject = null;

            if (file != null && properFilename in tableObjects)
                tableObject = tableObjects[properFilename]
            // else if (state.path[state.path.length - 1] == "Parameters")
            //     tableObject = tableObjects["Parameters"];

            if (file != null && properFilename in plotObjects)
                plotObject = plotObjects[properFilename]
    
            return {
                selected,
                properFilename,
                isFile: file != null,
                tableObject,
                plotObject
            };
        }
    },
    getProperFilenameOf: (state, getters) => {
        return (filename) => {
            if (filename == null)
                return null;

            const properFilenames = (Object.keys(state.outputFiles).length > 0) ? {
                ...state.inputFiles,
                ...state.outputFiles,
                ...state.outputFiles["Intermediate Files"].files,
                ...(state.outputFiles["Plottables"] || {files:[]}).files
            } : state.inputFiles;

            if (filename in properFilenames)
                return filename
            else if (filename == "target.0")
                return "target"
        
            const nnnFilename = filename.replaceAll(/\d{3}/g, "nnn");
            const isCFile = filename.slice(-2) == ".c";
            const isBswFile = filename.slice(-4) == ".bsw";
        
            if (nnnFilename in properFilenames)
                return nnnFilename;
            else if (isCFile)
                return "name.c"
            else if (isBswFile)
                return "name.bsw"
        
            return null;
        }
    }
}

export default {
    namespaced: true,
    state,
    getters,
    actions,
    mutations
}

// Helper functions
function getProperFilenameOf(filename, inputFiles, path) {
    if (path.length == 4)
        return path[2];

    if (filename in inputFiles)
        return filename
    else if (filename == "target.0")
        return "target"

    const nnnFilename = filename.replaceAll(/\d{3}/g, "nnn");
    const isCFile = filename.slice(-2) == ".c";
    const isBswFile = filename.slice(-4) == ".bsw";

    if (nnnFilename in inputFiles)
        return nnnFilename;
    else if (isCFile)
        return "name.c"
    else if (isBswFile)
        return "name.bsw"

    return null;
}

function dependenciesMatchesPath(dependencies, path) {
    return dependencies.some((dependency) => 
        dependency.every((label, i) => i >= path.length || path[i] == label))
}
