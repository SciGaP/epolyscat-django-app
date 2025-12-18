const state = {
    loading: {}
};

const actions = {
    
};

const mutations = {
    START(state, { key, message }) {
        state.loading = {...state.loading}; // This causes vue to update.
        state.loading[key] = state.loading[key] || {};
        state.loading[key][message] = {}; // state[key] is acting as a hashset the value message is asigned to is arbitrary
    },
    STOP(state, { key, message }) {
        state.loading = {...state.loading}; // This causes vue to update.
        state.loading[key] = state.loading[key] || {};
        delete state.loading[key][message];
    },
    CLEAR(state, { key }) {
        state.loading = {...state.loading}; // This causes vue to update.
        state.loading[key] = {};
    }
};

const getters = {
    getMessages(state) {
        return (key) => Object.keys(state.loading[key] || {});
    }
};

export default {
    namespaced: true,
    state,
    getters,
    actions,
    mutations,
};
