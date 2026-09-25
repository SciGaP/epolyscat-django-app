import { SettingsService } from "@/service/epolyscat-service";
let settingsRequest = null;

const state = {
  settings: {},
  loading: false,
  loaded: false,
  error: "",
  prefrences: {}
};

const actions = {
  async fetchSettings({ commit }) {
    if (settingsRequest) return settingsRequest;
    commit("settingsLoading");
    settingsRequest = (async () => {
      try {
        const settings = await SettingsService.all();
        commit("setSettings", settings);
        return true;
      } catch (error) {
        commit("settingsFailed", "Unable to load applications from Airavata. Please try again.");
        return false;
      } finally {
        settingsRequest = null;
      }
    })();
    return settingsRequest;
  },
};

const mutations = {
  setSettings(state, settings) {
    state.settings = settings;
    state.loading = false;
    state.loaded = true;
    state.error = "";
  },
  settingsLoading(state) {
    state.loading = true;
    state.loaded = false;
    state.settings = {};
    state.error = "";
  },
  settingsFailed(state, message) {
    state.loading = false;
    state.loaded = false;
    state.settings = {};
    state.error = message;
  },
  setPreference(state, {key, value}) {
       state.prefrences[key] = value;
  }
};

const getters = {
  epolyscatApplicationModuleId: (state) => {
    return state.settings && state.settings.EPOLYSCAT
      ? state.settings.EPOLYSCAT.EPOLYSCAT_APPLICATION_ID
      : null;
  },
  gaussian16ApplicationModuleId: (state) => {
    return state.settings && state.settings.EPOLYSCAT
      ? state.settings.EPOLYSCAT.GAUSSIAN16_APPLICATION_ID
      : null;
  },
  openmolcasApplicationModuleId: (state) => {
    return state.settings && state.settings.EPOLYSCAT
      ? state.settings.EPOLYSCAT.OPENMOLCAS_APPLICATION_ID
      : null;
  },
  getPreference: (state) => {
        return (key) => state.prefrences[key]
  }
};

export default {
  namespaced: true,
  state,
  getters,
  actions,
  mutations,
};
