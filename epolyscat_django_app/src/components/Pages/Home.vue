<template>
  <div class="portal-page">
    <main class="portal-page-content">
      <header class="portal-page-header"><h1 class="portal-page-title">ePolyScat</h1></header>
      <div class="portal-home-grid">
        <section aria-label="About ePolyScat">
          <p>ePolyScat is an electron–molecule scattering and photoionization suite for computing photoelectron spectra, differential cross sections, resonances, and angular distributions. This portal lets you create, submit, and analyze ePS runs on modern HPC resources.</p>
          <b-button variant="primary" @click="goToNewRun">New Run</b-button>
        </section>
        <section aria-labelledby="recent-runs-title">
          <h2 id="recent-runs-title">Recent Runs</h2>
          <LoadingOverlay name="runs">
            <div class="portal-recent-runs">
              <router-link v-for="run in displayedRuns" :key="run.id" :to="`/runs/${run.id}`" v-slot="{isExactActive, href, navigate}">
                <b-button variant="outline-primary" :class="{active: isExactActive}" :href="href" @click="navigate">{{ run.name }}</b-button>
              </router-link>
              <p v-if="runCount === 0" class="portal-empty-state">No recent runs</p>
            </div>
          </LoadingOverlay>
          <b-button v-if="runCount > 0" variant="outline-primary" to="/runs">View all runs</b-button>
        </section>
      </div>
    </main>
  </div>
</template>

<script>
//import ExperimentCard from "@/components/block/ExperimentCard";
import store from "@/store";
import LoadingOverlay from '../overlay/LoadingOverlay.vue';
import { eventBus } from '@/event-bus';

export default {
  components: { LoadingOverlay },
  name: 'Home',
  store: store,
  //components: { ExperimentCard },
  data() {
    return {
      //numberOfExperiments: 5,
      //experimentStatistics: {
      //  experimentCount: null,
      //  runCount: null
      //}
    }
  },
  computed: {
    //experiments() {
    //  return this.$store.getters["experiment/getExperiments"]({ pageSize: this.numberOfExperiments });
    //},
    displayedRuns() {
       let runs = this.$store.getters["run/getRuns"]();

       runs.sort((run1, run2) =>
          (new Date(run2.created)).valueOf() -
             (new Date(run1.created)).valueOf()
       );

       return runs.slice(0, 4);
    },
    runCount() {
                return this.$store.getters["run/getRuns"]().length;
    }
  },
  methods: {
    goToNewRun() {
      this.$router.push("/runs/new");
    },
    async refreshData() {
      try {
                    await this.$store.dispatch("run/fetchRuns", {});

                    this.$store.getters["run/getRuns"]();

                    this.$store.commit("loading/STOP", { key: "runs", message: "Fetching Runs" });
      } catch (error) {
                    eventBus.$emit("error", { name: `Error while trying to fetch the runs`, error });
      }
      //this.$store.dispatch("experiment/fetchExperiments", { pageSize: this.numberOfExperiments });
      //this.$store.dispatch("run/fetchRuns");

      //this.experimentStatistics = await ExperimentService.fetchExperimentStatistics();
    }
  },
  mounted() {
    this.refreshData();
    this.refreshDataInterval = setInterval(() => {
       this.refreshData();
    }, 45000);
  },
  beforeDestroy() {
     clearInterval(this.refreshDataInterval);
 }
}
</script>

<style scoped>
.cardList, .cardList > * {
        width: -webkit-fill-available;
        width: -moz-available;
        margin: 0;
}

.cardList > * {
        margin: 10px;
}

</style>
