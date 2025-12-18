import { SettingsService } from "@/service/epolyscat-service";

const state = {
  settings: {},
  prefrences: {}
};

const actions = {
  async fetchSettings({ commit }) {
    const settings = await SettingsService.all();

    commit("setSettings", settings);
  },
};

const mutations = {
  setSettings(state, settings) {
    state.settings = settings;
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
