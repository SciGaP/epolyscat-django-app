import {RunService} from "@/service/trecx-service";

const state = {
    runListMap: {},
    runListPaginationMap: {},
    runMap: {},

    viewableContentMap: {}
}

const actions = {
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

export default {
    namespaced: true,
    state,
    getters,
    actions,
    mutations
}
