import {ViewService} from "@/service/epolyscat-service";

const state = {
    //viewListMap: {},
    //viewListPaginationMap: {},
    viewMap: {}
}

const actions = {
    //async fetchViews({commit}) {
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
 
        //const views = await ViewService.fetchViews();
        //commit("SET_VIEW_MAP", { views });

        //return views;
    //},
    async deleteView({commit}, { viewId }) {
        await ViewService.deleteView(viewId);

        commit("REMOVE_VIEW", { viewId });
        commit("run/REMOVE_VIEW", { viewId }, { root: true });
    },


    /*
    async fetchViews({commit}) {
                             , {page = 1, pageSize = 1000, tutorials = false} = {
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
    */
    async fetchView({commit}, {viewId}) {
        const view = await ViewService.fetchView({viewId});

        commit("SET_VIEW_MAP", { views: [view] });

        return view;
    },
    async createView({commit}, {name, runIds}) {
        const view = await ViewService.createView({name, runIds});
        //const {viewId, owner, updated, created, deleted, type, activeRunCount, runCount, readonly} = view;
        commit("SET_VIEW_MAP", { views: [view] });

        //commit("SET_VIEW", {
        //    viewId, name, owner, updated, created, deleted, type, activeRunCount, runCount, readonly
        //});

        return view;
    },
    //async updateView({commit}, {viewId, name, runIds}) {
    async updateView({commit}, {viewId, name, runIds, override }) {
        //const view = await ViewService.updateView({viewId, name, runIds});
        const view = await ViewService.updateView({viewId, name, runIds, override });
        //const {owner, updated, created, deleted, type, activeRunCount, runCount, readonly} = view;

        commit("SET_VIEW_MAP", { views: [view] });

        //commit("SET_VIEW", {
        //    viewId, name, runIds, owner, updated, created, deleted, type, activeRunCount, runCount, readonly
        //});

        return view;
    },
    async insertIntoViews({commit, getters}, { viewIds, run }) {
        for (const viewId of viewIds) {
            let view = getters["getView"](viewId);

            if (!view) {
                view = await ViewService.fetchView(viewId);
                commit("SET_VIEW_MAP", { views: [view] });
            }

            commit("INSERT_RUN", { viewId, run });
        }
    }
}

const mutations = {
    UPDATE_VIEW(state, { viewId, viewData }) {
        Object.entries(viewData).forEach(([key, value]) => {
            state.viewMap[viewId][key] = value;
        });

        state.viewMap = {...state.viewMap};
    },
    SET_VIEW_MAP(state, { views }) {
        views.forEach(view => state.viewMap[view.id] = view);
        state.viewMap = {...state.viewMap};
    },
    REMOVE_VIEW(state, { viewId }) {
        delete state.viewMap[viewId];
        state.viewMap = {...state.viewMap};
    },
    INSERT_RUN(state, { viewId, run }) {
        // console.log(viewId, state.viewMap[viewId])
        if (!(viewId in state.viewMap)) {
          //print to log
          console.log('no viewId in vewMap')
        }

        state.viewMap[viewId].runs.push(run);
        state.viewMap[viewId].runCount = state.viewMap[viewId].runs.length;
    },
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

/*
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
*/
    getViews: (state) => {
        return () => Object.values(state.viewMap).filter(view => view.id != -1);
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
