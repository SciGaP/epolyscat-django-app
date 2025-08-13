import { SettingsService } from "@/service/trecx-service";

const state = {
  settings: {},
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
};

const getters = {
  trecxApplicationModuleId: (state) => {
    return state.settings && state.settings.TRECX
      ? state.settings.TRECX.TRECX_APPLICATION_ID
      : null;
  },
};

export default {
  namespaced: true,
  state,
  getters,
  actions,
  mutations,
};
