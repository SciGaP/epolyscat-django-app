import { plotObjects } from "@/fileData";
import { RunService, InputService } from "@/service/epolyscat-service";

const OUTPUTS = ["bound_tab", "H.DAT", "KMAT", "TMAT", "KMAT_keep", "TMAT_keep", "tr_nnn_nnn"];

const state = {
    //runListMap: {},
    //runListPaginationMap: {},
    runMap: {},

    //viewableContentMap: {}
}

const actions = {
    async fetchRuns({commit}) {
        const runs = await RunService.fetchRuns();

        commit("SET_RUN_MAP", { runs });

        return runs;
    },

    /*
    async fetchRuns({commit}, {experimentId = null, viewId = null, page = 1, pageSize = 1000} = {
        page: 1,
        pageSize: 1000
    }) {
        const queryString = JSON.stringify({experimentId, viewId, page, pageSize});

        const runsRes = await RunService.fetchAllRuns({experimentId, viewId, page, pageSize});
        const runs = runsRes.results;

        const runIds = runs.map((
            {
                runId, name, experimentId, created, updated, deleted, inpcDownloadUrl, groupResourceProfileId,
                computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory, status,
                resource, resourceShort, executions, canResubmit, canSubmit, inputTable
            }) => {

            commit("SET_RUN", {
                runId, name, experimentId, created, updated, deleted, inpcDownloadUrl, groupResourceProfileId,
                computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory, status,
                resource, resourceShort, executions, canResubmit, canSubmit, inputTable
            });

            return runId;
        });

        commit("SET_RUN_LIST", {
            queryString,
            runIds,
            pagination: {page, pageSize, total: runsRes.count}
        });
    },
    */
    async fetchRun({commit, getters}, { runId }) {
        let run;

        try {
            run = await RunService.fetchRun(runId);
        } catch (error) {
            if (getters["getRun"](runId) != undefined)
                run = getters["getRun"](runId)
            else
                throw new Error(`Error trying to find run with id: ${runId}, with error message: ${error}`)
        }

        commit("SET_RUN_MAP", { runs: [ run ] });

        return run;
    },
    /*
    async fetchRun({commit}, {runId = null} = {}) {
        const run = await RunService.fetchRun({runId});
        const {
            name, experimentId, created, updated, deleted, inpcDownloadUrl, groupResourceProfileId,
            computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory, status,
            resource, resourceShort, executions, canResubmit, canSubmit, inputTable
        } = run;

        commit("SET_RUN", {
            runId, name, experimentId, created, updated, deleted, inpcDownloadUrl, groupResourceProfileId,
            computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory, status,
            resource, resourceShort, executions, canResubmit, canSubmit, inputTable
        });
    }, 
    */
    async loadInputs({commit, getters, dispatch}, { runId }) {
        let run;

        try {
            run = getters["getRun"](runId) || await dispatch("fetchRun");
            run.inputs.length; // ensures that run is not undefined
        } catch (error) {
            throw new Error(`Error trying to find run with id: ${runId}, with error message: ${error}`);
        }

        let path = [];

        for (const input of run.inputs) {
            if (input.type == "files")
                for (const file of input.files) {
                    file.unchanged = true;

                    commit("input/ADD_TO_INPUT_FILE", {
                        file,
                        inputFileName: input.name
                    }, { root: true })
                }
            else if (input.type == "parameter")
                commit("input/SET_PARAMETER", {
                    name: input.name,
                    value: input.value
                }, { root: true });
            else if (input.name == "Calculation_Type")
                path[0] = input.value;
            else if (["EPOLYSCAT_Application_Module", "Application_Utility", "Application_Workflow"].indexOf(input.name) != -1)
                path[1] = input.value
        }

        commit("input/SET_PATH", { path }, { root: true });

        return run;
    },
    async tryLoadOutputs({commit, getters, dispatch, rootGetters}, { runId }) {
        let run;

        const path = rootGetters["input/getPath"];
        console.log("Path being gotten: ", path);

        try {
            run = getters["getRun"](runId) || await dispatch("fetchRun");
            run.inputs.length; // ensures that run is not undefined
        } catch (error) {
            throw new Error(`Error trying to find run with id: ${runId}, with error message: ${error}`);
        }

        const outputFileObjs = await InputService.fetchOutputs(runId);

        if (outputFileObjs.length == 0)
            return run;

        const inputFiles = run.inputs
            .filter(input => input.type == "files")
            .reduce((acc, input) => [...acc, ...input.files], []);
        let outputFiles = {};

        for (const inputFile of inputFiles) {
            if (inputFile != null && typeof inputFile == "object" && !("dataProductURI" in inputFile))
                inputFile.dataProductURI = inputFile["data-product-uri"]

            const category = "Input Files";
            const secondaryCategory =
                (/\.c$/.test(inputFile.name)) ? "name.c" :
                (/\.bsw$/.test(inputFile.name)) ? "name.bsw" :
                (/\d{3}/.test(inputFile.name)) ? inputFile.name.replaceAll(/\d{3}/g, "nnn") :
                inputFile.name;

            outputFiles[category] = outputFiles[category] || {
                // isMultiFileInput: false,
                files: {},
                isValid: () => true,
                dependencies: [path.slice(0,2)]
            };

            outputFiles[category].files[secondaryCategory] = outputFiles[category].files[secondaryCategory] || {
                // isMultiFileInput: false,
                files: [],
                isValid: () => true,
                dependencies: [[...path.slice(0,2), category]]
            };

            outputFiles[category].files[secondaryCategory].files.push(inputFile);
            outputFiles[category].files[secondaryCategory].files.sort((a,b) => a.name.localeCompare(b.name));
        }

        for (const outputFileObj of outputFileObjs) {
            if (outputFileObj != null && typeof outputFileObj == "object" && !("dataProductURI" in outputFileObj))
                outputFileObj.dataProductURI = outputFileObj["dataProductURI"]

            const isInput = inputFiles.some(file => file.name == outputFileObj.name);

            const isJobFile = /slurm|A\d{5,}/.test(outputFileObj.name);
            const isPlottable = outputFileObj.name in plotObjects ||
                outputFileObj.name.replaceAll(/\d{3}/g, "nnn") in plotObjects;

            if (OUTPUTS.indexOf(outputFileObj.name) >= 0/* && !isInput*/) { // Assumes that there is only ever one output file
                outputFiles[outputFileObj.name] = {
                    // isMultiFileInput: false,
                    files: [outputFileObj],
                    isValid: () => true,
                    dependencies: [path.slice(0,2)]
                };
            } else if (isJobFile /*&& !isInput*/) {
                outputFiles["Job Info"] = outputFiles["Job Info"] || {
                    // isMultiFileInput: false,
                    files: [],
                    isValid: () => true,
                    dependencies: [path.slice(0,2)]
                };

                outputFiles["Job Info"].files.push(outputFileObj);
                outputFiles["Job Info"].files.sort((a,b) => a.name.localeCompare(b.name));
            } else {
                const category =
                    (isPlottable) ? "Plottables" :
                    (outputFileObj.name.includes("log")) ? "Logs" :
                    "Intermediate Files";

                const secondaryCategory =
                    (isInput && /\.c$/.test(outputFileObj.name)) ? "name.c" :
                    (isInput && /\.bsw$/.test(outputFileObj.name)) ? "name.bsw" :
                    (/\d{3}/.test(outputFileObj.name)) ? outputFileObj.name.replaceAll(/\d{3}/g, "nnn") :
                    outputFileObj.name;

                outputFiles[category] = outputFiles[category] || {
                    // isMultiFileInput: false,
                    files: {},
                    isValid: () => true,
                    dependencies: [path.slice(0,2)]
                };

                outputFiles[category].files[secondaryCategory] = outputFiles[category].files[secondaryCategory] || {
                    // isMultiFileInput: false,
                    files: [],
                    isValid: () => true,
                    dependencies: [[...path.slice(0,2), category]]
                };

                outputFiles[category].files[secondaryCategory].files.push(outputFileObj);
                outputFiles[category].files[secondaryCategory].files.sort((a,b) => a.name.localeCompare(b.name));
            }
        }
             
        commit("input/SET_INPUTS", { outputFiles }, { root: true });

        return run;
    },
    async createRun({dispatch, commit, rootGetters}, {
        name, groupResourceProfileId, computeResourceId, coreCount,
        totalPhysicalMemory, nodeCount, wallTimeLimit, queueName, viewIds, description
    }) {
        const inputs = await rootGetters["input/getPreparedInputs"]({ prepareForCreation: true });

        const data = await RunService.createRun({
            name, inputs, groupResourceProfileId, computeResourceId,
            coreCount, nodeCount, wallTimeLimit, queueName, totalPhysicalMemory, viewIds, description
        }, false);

        await dispatch("view/insertIntoViews", { viewIds, run: data }, { root: true });

        return data;
    },
    async updateRun({rootGetters}, {
        name, id, groupResourceProfileId, computeResourceId, coreCount,
        totalPhysicalMemory, nodeCount, wallTimeLimit, queueName, description
    }) {
        const inputs = await rootGetters["input/getPreparedInputs"]({ prepareForCreation: false });

        const data = await RunService.updateRun({
            name, runId: id, inputs, groupResourceProfileId, computeResourceId,
            coreCount, nodeCount, wallTimeLimit, queueName, totalPhysicalMemory, description
        }, false);

        return data;
    },
    async submitRun({commit}, { runId }) {
        const run = await RunService.submitRun({ runId });

        run.status = "EXECUTING";
        run.jobStatus = "NO STATUS";
        run.displayStatus = "EXECUTING";

        commit("SET_RUN_MAP", { runs: [run] });

        return run;
    },
    async deleteRun({commit}, { runId, deleteAssociated }) {
        await RunService.deleteRun({ runId, deleteAssociated });

        commit("REMOVE_RUN", { runId });
    },
    async cloneRun({commit, rootGetters}, {
        name, groupResourceProfileId, computeResourceId, coreCount,
        totalPhysicalMemory, nodeCount, wallTimeLimit, queueName, inputs, description
    }) {
        if (inputs == null)
            inputs = await rootGetters["input/getPreparedInputs"]({ prepareForCreation: true });

        name = `Clone of ` + name;

        const run = await RunService.createRun({
            name, inputs, groupResourceProfileId, computeResourceId,
            coreCount, nodeCount, wallTimeLimit, queueName, totalPhysicalMemory, description
        }, false);

        // console.log("NEW RUN ID" + run.id)

        commit("SET_RUN_MAP", { runs: [ run ] });

        return run;
    },
    async fetchStatus({}, {runId}) {
        return await RunService.fetchStatus({runId});
    },
    async changeNotificationSettings({commit}, {runId, isEmailNotificationOn}) {
        const run = RunService.changeNotificationSettings({ runId, isEmailNotificationOn });

        commit("SET_RUN_MAP", { runs: [ run ] });

        return run;
    }
}
/*
    async createRun({commit}, {root, experimentId, directedit, groupResourceProfileId, computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory}) {
        const {name, runId} = await RunService.createRun({
            root, experimentId, directedit, groupResourceProfileId, computeResourceId, queueName, coreCount, nodeCount,
            wallTimeLimit, totalPhysicalMemory
        });

        commit("SET_RUN", {
            runId, name, experimentId, directedit, groupResourceProfileId, computeResourceId, queueName, coreCount,
            nodeCount, wallTimeLimit, totalPhysicalMemory
        });
    },
    setRun({commit}, {
        runId, name, experimentId, created, updated, deleted, inpcDownloadUrl, groupResourceProfileId,
        computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory, status,
        resource, resourceShort, executions, canResubmit, canSubmit, inputTable
    }) {
        commit("SET_RUN", {
            runId, name, experimentId, created, updated, deleted, inpcDownloadUrl, groupResourceProfileId,
            computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory, status,
            resource, resourceShort, executions, canResubmit, canSubmit, inputTable
        });
    },
    async fetchViewableContent({commit}, {runId, filename, inpcDownloadUrl}) {
        const content = await RunService.fetchViewableContent({runId, filename, inpcDownloadUrl});
        commit("SET_VIEWABLE_CONTENT", {runId, filename, inpcDownloadUrl, content});
    }
}

const mutations = {
    SET_RUN_LIST(state, {queryString, runIds, pagination: {page, pageSize, total}}) {
        state.runListMap = {
            ...state.runListMap,
            [queryString]: runIds
        };
        state.runListPaginationMap = {
            ...state.runListPaginationMap,
            [queryString]: {page, pageSize, total}
        };
    },
    SET_RUN(state, {
        runId, name, experimentId, created, updated, deleted, inpcDownloadUrl, groupResourceProfileId,
        computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory, status,
        resource, resourceShort, executions, canResubmit, canSubmit, inputTable
    }) {
        state.runMap = {
            ...state.runMap,
            [runId]: {
                runId, name, experimentId, created, updated, deleted, inpcDownloadUrl, groupResourceProfileId,
                computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory, status,
                resource, resourceShort, executions, canResubmit, canSubmit, inputTable
            }
        }
    },
    SET_VIEWABLE_CONTENT(state, {runId, filename, inpcDownloadUrl, content}) {
        state.viewableContentMap = {
            ...state.viewableContentMap,
            [`${runId}-${filename}-${inpcDownloadUrl}`]: content
        }
    }
}


const getters = {
    getRuns: (state, getters) => {
        return ({experimentId = null, viewId = null, page = 1, pageSize = 1000} = {page: 1, pageSize: 1000}) => {
            const queryString = JSON.stringify({experimentId, viewId, page, pageSize});

            const runIds = state.runListMap[queryString];
            if (runIds) {
                return runIds.map(runId => getters.getRun({runId}));
            } else {
                return null;
            }
        }
    },
    getRunsPagination: (state) => {
        return ({experimentId = null, viewId = null, page = 1, pageSize = 1000} = {page: 1, pageSize: 1000}) => {
            const queryString = JSON.stringify({experimentId, viewId, page, pageSize});
            const runListPagination = state.runListPaginationMap[queryString];
            if (runListPagination) {
                return runListPagination;
            } else {
                return null;
            }
        }
    },
    getRun: (state) => {
        return ({runId}) => {
            if (state.runMap[runId]) {
                return state.runMap[runId];
            } else {
                return null;
            }
        }
    },
    getViewableContent: (state) => {
        return ({runId, filename, inpcDownloadUrl}) => {
            if (state.viewableContentMap[`${runId}-${filename}-${inpcDownloadUrl}`]) {
                return state.viewableContentMap[`${runId}-${filename}-${inpcDownloadUrl}`];
            } else {
                return null;
            }
        }
    }
}
*/
const mutations = {
    UPDATE_RUN(state, { runId, runData }) {
        Object.entries(runData).forEach(([key, value]) => {
            state.runMap[runId][key] = value;
        });

        state.runMap = {...state.runMap}; // Ensures that vuex knows that the value was updated
    },
    SET_RUN_MAP(state, { runs }) {
        runs.forEach(run => state.runMap[run.id] = run);
        state.runMap = {...state.runMap};
        // console.log("STORE: runs has been modified", state.runMap);
    },
    REMOVE_RUN(state, { runId }) {
        delete state.runMap[runId];
        state.runMap = {...state.runMap};
    },
    REMOVE_VIEW(state, { viewId }) {
        Object.keys(state.runMap).forEach(runId => {
            state.runMap[runId].viewIds = state.runMap[runId].viewIds.filter(otherViewId => otherViewId == viewId);
        })
    }
}

const getters = {
    getRuns: (state) => {
        return () => Object.values(state.runMap).filter(run => !run.isTutorial);
    },
    getRun: (state) => {
        return (runId) => state.runMap[runId];
    },
    getFirstNRuns: (state) => {
        return (n) => {
            let runs = Object.values(state.runMap);
            runs.sort((a, b) => a.id - b.id);
            return runs.slice(0, n);
        }
    },
    getStrippedRun: (state) => {
        return (run) => {
            return {
                name: run.name,
                id: run.id,
                isEmailNotificationOn: run.isEmailNotificationOn,
                description: run.description,
                computeResourceId: run.computeResourceId,
                queueName: run.queueName,
                coreCount: run.coreCount,
                nodeCount: run.nodeCount,
                wallTimeLimit: run.wallTimeLimit,
                totalPhysicalMemory: run.totalPhysicalMemory,
                inputs: run.inputs
            };
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
