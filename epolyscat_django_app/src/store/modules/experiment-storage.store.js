import {ExperimentService} from "@/service/trecx-service";

const state = {
    experimentListMap: {},
    experimentListPaginationMap: {},
    experimentMap: {}
}

const actions = {
    async fetchExperiments({commit}, {page = 1, pageSize = 1000} = {page: 1, pageSize: 1000}) {
        const queryString = JSON.stringify({page, pageSize});

        const experimentsRes = await ExperimentService.fetchAllExperiments({page, pageSize});
        const experiments = experimentsRes.results;

        const experimentIds = experiments.map(({experimentId, name, root, description, activeRunCount, runCount, owner, updated, created, deleted, airavataProjectId}) => {
            commit("SET_EXPERIMENT", {
                experimentId, name, root, description, activeRunCount, runCount, owner,
                updated, created, deleted, airavataProjectId
            });

            return experimentId;
        });

        commit("SET_EXPERIMENT_LIST", {
            queryString,
            experimentIds,
            pagination: {page, pageSize, total: experimentsRes.count}
        });
    },
    async fetchExperiment({commit}, {experimentId}) {
        const experiment = await ExperimentService.fetchExperiment({experimentId});
        const {name, root, description, activeRunCount, runCount, owner, updated, created, deleted, airavataProjectId} = experiment;
        commit("SET_EXPERIMENT", {
            experimentId, name, root, description, activeRunCount, runCount, owner, updated, created,
            deleted, airavataProjectId
        });
    },
    async createExperiment({commit}, {name, description}) {
        const experiment = await ExperimentService.createExperiment({name, description});
        const {experimentId, root, activeRunCount, runCount, owner, updated, created, deleted, airavataProjectId} = experiment;
        commit("SET_EXPERIMENT", {
            experimentId, name, root, description, activeRunCount, runCount, owner,
            updated, created, deleted, airavataProjectId
        });

        return experiment;
    }
}

const mutations = {
    SET_EXPERIMENT_LIST(state, {queryString, experimentIds, pagination: {page, pageSize, total}}) {
        state.experimentListMap = {
            ...state.experimentListMap,
            [queryString]: experimentIds
        }
        state.experimentListPaginationMap = {
            ...state.experimentListPaginationMap,
            [queryString]: {page, pageSize, total}
        };
    },
    SET_EXPERIMENT(state, {experimentId, id, name, root, description, activeRunCount, runCount, owner, updated, created, deleted, airavataProjectId}) {
        state.experimentMap = {
            ...state.experimentMap,
            [experimentId]: {
                experimentId, id, name, root, description, activeRunCount, runCount, owner, updated, created,
                deleted, airavataProjectId
            }
        }
    }
}


const getters = {
    getExperiments: (state, getters) => {
        return ({page = 1, pageSize = 1000} = {page: 1, pageSize: 1000}) => {
            const queryString = JSON.stringify({page, pageSize});
            const experimentIds = state.experimentListMap[queryString];
            if (experimentIds) {
                return experimentIds.map(experimentId => getters.getExperiment({experimentId}));
            } else {
                return null;
            }
        }
    },
    getExperimentsPagination: (state) => {
        return ({page = 1, pageSize = 1000} = {page: 1, pageSize: 1000}) => {
            const queryString = JSON.stringify({page, pageSize});
            const experimentListPagination = state.experimentListPaginationMap[queryString];
            if (experimentListPagination) {
                return experimentListPagination;
            } else {
                return null;
            }
        }
    },
    getExperiment: (state) => {
        return ({experimentId}) => {
            if (state.experimentMap[experimentId]) {
                return state.experimentMap[experimentId];
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