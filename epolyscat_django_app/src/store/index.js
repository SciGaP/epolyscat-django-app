import Vue from 'vue';
import Vuex from 'vuex';
import createLogger from 'vuex/dist/logger';

import experimentStore from './modules/experiment-storage.store';
import runStore from './modules/run-storage.store';
import viewStore from './modules/view-storage.store';
import settingsStore from './modules/settings.store';

Vue.use(Vuex);

const debug = true;

export default new Vuex.Store({
    modules: {
        "experiment": experimentStore,
        "run": runStore,
        "view": viewStore,
        settings: settingsStore,
    },
    strict: debug,
    plugins: debug ? [createLogger()]: [],
});
