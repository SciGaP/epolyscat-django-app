<template>
  <div class="w-100 h-100 overflow-auto">
    <div class="w-100 d-flex p-5">
      <div class="flex-fill">
        <div class="text-primary" style="font-size:56px;">ePolyScat Manual</div>
        <div>
          ePolyScat is an electron–molecule scattering and photoionization
          suite for computing photoelectron spectra, differential cross sections,
          resonances, and angular distributions. This portal lets you create,
          submit, and analyze ePS runs on modern HPC resources.
        </div>
      </div>
      <div class="text-primary text-center p-5 m-2 stats-card">
        <div style="font-size:56px;height:85px;">
          <template v-if="experimentStatistics">{{ experimentStatistics.experimentCount }}</template>
          <template v-else>0</template>
        </div>
        <div>Tutorials</div>
      </div>
      <div class="text-primary text-center p-5 m-2 stats-card">
        <div style="font-size:56px;height:85px;">
          <template v-if="experimentStatistics">{{ experimentStatistics.runCount }}</template>
          <template v-else>0</template>
        </div>
        <div>Runs</div>
      </div>
    </div>
    <div class="exp-list p-3">
      <div class="grid-container">
        <router-link to="/create-experiment" v-slot="{ href, route, navigate, isActive,isExactActive }">
          <b-link :href="href" @click="navigate">
            <div class="exp-card create-exp-card p-4 text-primary text-center">
              <div style="font-size: 50px;">
                <b-icon icon="plus-circle"/>
              </div>
              <div style="font-size: 20px;">Create a new experiment</div>
            </div>
          </b-link>
        </router-link>

        <!-- experiments loop -->
        <div v-for="experiment in experiments" :key="experiment.experimentId" class="exp-card">
          <ExperimentCard :experimentId="experiment.experimentId"/>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import ExperimentCard from "@/components/block/ExperimentCard";
import store from "@/store";
import { ExperimentService } from "@/service/epolyscat-service";

export default {
  name: 'Home',
  store: store,
  components: { ExperimentCard },
  data() {
    return {
      numberOfExperiments: 5,
      experimentStatistics: {
        experimentCount: null,
        runCount: null
      }
    }
  },
  computed: {
    experiments() {
      return this.$store.getters["experiment/getExperiments"]({ pageSize: this.numberOfExperiments });
    }
  },
  methods: {
    async refreshData() {
      this.$store.dispatch("experiment/fetchExperiments", { pageSize: this.numberOfExperiments });
      this.$store.dispatch("run/fetchRuns");

      this.experimentStatistics = await ExperimentService.fetchExperimentStatistics();
    }
  },
  mounted() {
    this.refreshData();
  }
}
</script>

<style scoped>
.stats-card {
  background: #2265970D;
}

.exp-card {
  width: 300px !important;
  height: 180px !important;
}

.create-exp-card {
  border: 2px dashed #1C6296;
  box-sizing: border-box;
  border-radius: 5px;
  font-weight: bold;
  line-height: 70px;
}

.grid-container {
  column-gap: 20px !important;
  row-gap: 40px !important;
  display: grid;
  grid-auto-columns: auto;
  grid-template-columns: repeat( auto-fit, 300px);
}
</style>
