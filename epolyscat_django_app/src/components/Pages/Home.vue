<template>
  <!-- div class="w-100 h-100 overflow-auto" -->
  <div class="w-100 h-100 flex-fill d-flex flex-row p-5">

    <!-- div class="w-100 d-flex p-5" -->
      <!-- div class="flex-fill" -->
      <div style="padding: 0 10px; margin: 0 10px; width: 200%; height: 100%; overflow-y: scroll;">
        <div class="text-primary" style="font-size:56px;">ePolyScat</div>
        <hr>
        <p>
          ePolyScat is an electron–molecule scattering and photoionization
          suite for computing photoelectron spectra, differential cross sections,
          resonances, and angular distributions. This portal lets you create,
          submit, and analyze ePS runs on modern HPC resources.
        </p>
      </div>
      <div class="cardList m-5">
            <h3>Recent Runs</h3>
            <LoadingOverlay name="runs" class="cardList p-3" style="border-top: 1px solid #0002; border-bottom: 1px solid #0002">
                <router-link v-for="run in displayedRuns" :key="run" :to="`/runs/${run.id}`" v-slot="{isExactActive, href, navigate}">
                    <b-button variant="outline-primary" :class="{active: isExactActive}" :href="href" @click="navigate">
                        {{ run.name }}
                    </b-button>
                </router-link>
                <div v-if="runCount == 0" class="d-flex align-items-center justify-content-center" style="height: 200px">
                    <span>No recent runs</span>
                </div>
            </LoadingOverlay>
            <router-link :to="`runs/new`" v-slot="{isExactActive, href, navigate}">
                <b-button variant="outline-secondary" :class="{active: isExactActive}" :href="href" @click="navigate">
                    Create new run
                </b-button>
            </router-link>
            <router-link :to="`runs`" v-slot="{isExactActive, href, navigate}" v-if="runCount > 4">
                <b-button variant="outline-secondary" :class="{active: isExactActive}" :href="href" @click="navigate">
                    View all Runs
                </b-button>
            </router-link>
          <!-- router-link :to="`runs/new`" v-slot="{isExactActive, href, navigate}">
              <b-button variant="outline-secondary" :class="{active: isExactActive}" :href="href" @click="navigate">
                    Create new run
              </b-button>
          </router-link -->
      </div>
      
      <!-- div class="text-primary text-center p-5 m-2 stats-card">
        <div style="font-size:56px;height:85px;">
          <template v-if="experimentStatistics">{{ experimentStatistics.experimentCount }}</template>
          <template v-else>0</template>
        </div>
        <div>Experiments

        </div>
      </div>
      <div class="text-primary text-center p-5 m-2 stats-card">
        <div style="font-size:56px;height:85px;">
          <template v-if="experimentStatistics">{{ experimentStatistics.runCount }}</template>
          <template v-else>0</template>
        </div>
        <div>Runs

        </div>
      </div -->
    <!-- /div -->
    <!-- div class="exp-list p-3" -->
      <!-- div class="grid-container" -->
        <!-- router-link to="/create-experiment" v-slot="{ href, route, navigate, isActive,isExactActive }">
          <b-link :href="href" @click="navigate">
            <div class="exp-card create-exp-card p-4 text-primary text-center">
              <div style="font-size: 50px;">
                <b-icon icon="plus-circle"/>
              </div>
              <div style="font-size: 20px;">Create a new experiment</div>
            </div>
          </b-link>
        </router-link -->

        <!-- experiments loop -->
        <!-- div v-for="experiment in experiments" :key="experiment.experimentId" class="exp-card">
          <ExperimentCard :experimentId="experiment.experimentId"/>
        </div -->
      <!-- /div>
    </div -->
  </div>
</template>

<script>
//import ExperimentCard from "@/components/block/ExperimentCard";
import store from "@/store";
import LoadingOverlay from '../overlay/LoadingOverlay.vue';
import { eventBus } from '@/event-bus';
import { ExperimentService } from "@/service/epolyscat-service";

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
