import {ViewService} from "@/service/epolyscat-service";

const pendingViewRequests = {};

const state = {
    // Must be declared here so they are reactive and never `undefined`.
    // getViewsPagination reads viewListPaginationMap during the very first
    // render (before fetchViews resolves); leaving it undefined throws a
    // TypeError that aborts the Views page render (list shows only on the
    // second visit, once a prior fetch has lazily created the property).
    viewListMap: {},
    viewListPaginationMap: {},
    viewMap: {}
}

const actions = {
    //async fetchViews({commit}) {
    fetchViews({commit}, {page = 1, pageSize = 1000, tutorials = false} = {
        page: 1, pageSize: 1000, tutorials: false
    }) {
        const queryString = JSON.stringify({page, pageSize, tutorials});
        if (pendingViewRequests[queryString])
            return pendingViewRequests[queryString];

        const request = ViewService.fetchAllViews({page, pageSize, tutorials})
            .then(viewsRes => {
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
            });

        pendingViewRequests[queryString] = request.then(
            result => {
                delete pendingViewRequests[queryString];
                return result;
            },
            error => {
                delete pendingViewRequests[queryString];
                throw error;
            }
        );

        return pendingViewRequests[queryString];
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
            let view = getters["getView"]({viewId});

            if (!view) {
                view = await ViewService.fetchView({viewId});
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
        if (!(viewId in state.viewMap))
            return;
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
                // Return a safe default (not null) so the template's
                // `viewsPagination.total` never dereferences null before
                // fetchViews has populated the pagination map.
                return {page, pageSize, total: 0};
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
