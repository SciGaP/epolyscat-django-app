<template>
  <div class="exp-card w-100 h-100 p-4">
    <div class="exp-content w-100">
      <div class="w-100 overflow-auto" style="font-size: 20px;">{{ experiment.name }}</div>
      <div style="color: #989797;font-size: 12px;">Last edited on {{ experiment.updated }}</div>
      <div style="font-size: 12px;">
        <span style="color: #989797;">Runs: </span>
        <span>{{ experiment.activeRunCount }}/{{ experiment.runCount }} active runs</span>
      </div>
    </div>
    <router-link :to="`/runs/?experimentId=${experiment.experimentId}`"
                 v-slot="{ href, route, navigate, isActive,isExactActive }">
      <button class="exp-button w-100 p-2 text-primary mt-4 ml-1 mr-1" @click="navigate">View experiment</button>
    </router-link>
  </div>
</template>

<script>
import store from "@/store";

export default {
  name: 'ExperimentCard',
  store: store,
  props: ["experimentId"],
  computed: {
    experiment() {
      return this.$store.getters["experiment/getExperiment"]({experimentId: this.experimentId});
    }
  }
}
</script>

<style scoped>
.exp-card {
  background: #FFFFFF;
  border: 0.5px solid #D8D8D8;
  box-sizing: border-box;
  border-radius: 5px;
}

span {
  padding-bottom: 5px !important;
}

.exp-button-div {
  width: 100%;
  margin-left: -5px;
  margin-top: 30px;
}

.exp-button {
  /*width: 208px !important;*/
  /*border-radius: 5px;*/
  font-weight: 600;
  /*padding-top: 14px;*/
  /*padding-bottom: 14px;*/
  border: none;
  background: #1C62960D;
}
</style>