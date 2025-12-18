import Vue from 'vue';
import Vuex from 'vuex';
import createLogger from 'vuex/dist/logger';

import experimentStore from './modules/experiment-storage.store';
import RunStore from './modules/run-storage.store';
import ViewStore from './modules/view-storage.store';
import SettingsStore from './modules/settings.store';
import InputStore from './modules/input-storage.store';
import LoadingStore from './modules/loading.store';
import PlotParameters from './modules/plotParameters-storage.store';


Vue.use(Vuex);

const debug = true;

export default new Vuex.Store({
    modules: {
        "run": RunStore,
        "view": ViewStore,
        "input": InputStore,
        "settings": SettingsStore,
        "loading": LoadingStore,
        "plotParameters": PlotParameters

        //"experiment": experimentStore,
        //"run": runStore,
        //"view": viewStore,
        //settings: settingsStore,
    },
    strict: debug,
    plugins: debug ? [createLogger()]: [],
});
