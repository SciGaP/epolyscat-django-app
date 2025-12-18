import { PlotService } from "@/service/epolyscat-service";

const state = {
    plotParameters: []
}

const actions = {
    async fetchPlotParameters({commit}) {
        const plotParameters = await PlotService.getPlotParameters();

        commit("UPDATE_PLOT_PARAMETERS", { plotParameters });

        return plotParameters;
    },
}

const mutations = {
    UPDATE_PLOT_PARAMETERS(state, { plotParameters }) {
        state.plotParameters = plotParameters;
    },
}

const getters = {
    getPlotParameters: (state) => {
        return state.plotParameters;
    },
}

export default {
    namespaced: true,
    state,
    getters,
    actions,
    mutations
}
