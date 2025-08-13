import {ViewService} from "@/service/trecx-service";

const state = {
    viewListMap: {},
    viewListPaginationMap: {},
    viewMap: {}
}

const actions = {
    async fetchViews({commit}, {page = 1, pageSize = 1000, tutorials = false} = {
        page: 1, pageSize: 1000, tutorials: false
    }) {
        const queryString = JSON.stringify({page, pageSize, tutorials});

        const viewsRes = await ViewService.fetchAllViews({page, pageSize, tutorials});
        const views = viewsRes.results;

        const viewIds = views.map(({viewId, name, owner, updated, created, deleted, type, activeRunCount, runCount, readonly}) => {
            commit("SET_VIEW", {
                viewId, name, owner, updated, created, deleted, type, activeRunCount, runCount, readonly
            });

            return viewId;
        });

        commit("SET_VIEW_LIST", {
            queryString,
            viewIds,
            pagination: {page, pageSize, total: viewsRes.count}
        });
    },
    async fetchView({commit}, {viewId}) {
        const view = await ViewService.fetchView({viewId});
        const {name, owner, updated, created, deleted, type, activeRunCount, runCount, readonly} = view;
        commit("SET_VIEW", {
            viewId, name, owner, updated, created, deleted, type, activeRunCount, runCount, readonly
        });
    },
    async createView({commit}, {name, runIds}) {
        const view = await ViewService.createView({name, runIds});
        const {viewId, owner, updated, created, deleted, type, activeRunCount, runCount, readonly} = view;
        commit("SET_VIEW", {
            viewId, name, owner, updated, created, deleted, type, activeRunCount, runCount, readonly
        });

        return view;
    },
    async updateView({commit}, {viewId, name, runIds}) {
        const view = await ViewService.updateView({viewId, name, runIds});
        const {owner, updated, created, deleted, type, activeRunCount, runCount, readonly} = view;
        commit("SET_VIEW", {
            viewId, name, runIds, owner, updated, created, deleted, type, activeRunCount, runCount, readonly
        });

        return view;
    }
}

const mutations = {
    SET_VIEW_LIST(state, {queryString, viewIds, pagination: {page, pageSize, total}}) {
        state.viewListMap = {
            ...state.viewListMap,
            [queryString]: viewIds
        }
        state.viewListPaginationMap = {
            ...state.viewListPaginationMap,
            [queryString]: {page, pageSize, total}
        };
    },
    SET_VIEW(state, {viewId, name, owner, updated, created, deleted, type, activeRunCount, runCount, readonly}) {
        state.viewMap = {
            ...state.viewMap,
            [viewId]: {viewId, name, owner, updated, created, deleted, type, activeRunCount, runCount, readonly}
        }
    }
}


const getters = {
    getViews: (state, getters) => {
        return ({page = 1, pageSize = 1000, tutorials = false} = {page: 1, pageSize: 1000, tutorials: false}) => {
            const queryString = JSON.stringify({page, pageSize, tutorials});
            const viewIds = state.viewListMap[queryString];
            if (viewIds) {
                return viewIds.map(viewId => getters.getView({viewId}));
            } else {
                return null;
            }
        }
    },
    getViewsPagination: (state) => {
        return ({page = 1, pageSize = 1000, tutorials = false} = {page: 1, pageSize: 1000, tutorials: false}) => {
            const queryString = JSON.stringify({page, pageSize, tutorials});
            const viewListPagination = state.viewListPaginationMap[queryString];
            if (viewListPagination) {
                return viewListPagination;
            } else {
                return null;
            }
        }
    },
    getView: (state) => {
        return ({viewId}) => {
            if (state.viewMap[viewId]) {
                return state.viewMap[viewId];
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