<template>
  <div class="w-100 h-100 bg-light p-2">
    <div class="h-100 bg-white overflow-auto p-3" style="flex-grow: 1">
      <div class="w-100">
        <b-breadcrumb>
          <template v-if="experimentId">
            <b-breadcrumb-item to="/experiments">Experiments</b-breadcrumb-item>
            <b-breadcrumb-item v-if="experiment" :to="`/runs?experimentId=${experiment.experimentId}`">
              {{ experiment.name }}
            </b-breadcrumb-item>
          </template>
          <template v-else-if="viewId">
            <b-breadcrumb-item to="/Views">Views</b-breadcrumb-item>
            <b-breadcrumb-item v-if="view" :to="`/runs?viewId=${view.viewId}`">
              {{ view.name }}
            </b-breadcrumb-item>
          </template>
          <b-breadcrumb-item v-else to="/runs">Runs</b-breadcrumb-item>

          <template v-if="runId">
            <b-breadcrumb-item :to="runLink" v-if="run">
              {{ run.name }}
            </b-breadcrumb-item>
          </template>

          <template v-if="filename">
            <b-breadcrumb-item :to="viewableLink" v-if="filename">
              {{ filename }}
            </b-breadcrumb-item>
          </template>
        </b-breadcrumb>
        <Errors :errors="errors"/>
      </div>
      <div class="w-100 overflow-auto">
        <RunViewableEditor :run-id="runId" :run="run" :filename="filename" :read-only="true"/>
      </div>
    </div>
  </div>
</template>

<script>
import store from "@/store";
import Errors from "@/components/block/Errors";
import RunViewableEditor from "@/components/block/RunViewableEditor";

export default {
  name: 'CreateRun',
  components: {RunViewableEditor, Errors},
  store: store,
  data() {
    return {
      processing: false,

      errors: []
    };
  },
  computed: {
    viewableLink() {
      let _link = `/runs/${this.runId}/viewables/${this.filename}?`;

      if (this.experimentId) {
        _link += `experimentId=${this.experimentId}&`;
      }

      if (this.viewId) {
        _link += `viewId=${this.viewId}&`;
      }

      return _link;
    },
    runLink() {
      let _link = "/runs/";
      if (this.runId) {
        _link += `${this.runId}?`;
      }

      if (this.experimentId) {
        _link += `experimentId=${this.experimentId}&`;
      }

      if (this.viewId) {
        _link += `viewId=${this.viewId}&`;
      }

      return _link;
    },
    experimentId() {
      return this.$route.query.experimentId;
    },
    viewId() {
      return this.$route.query.viewId;
    },
    runId() {
      return this.$route.params.runId;
    },
    filename() {
      return this.$route.params.filename;
    },
    run() {
      return this.$store.getters["run/getRun"]({
        runId: this.runId
      });
    },
    experiment() {
      return this.$store.getters["experiment/getExperiment"]({
        experimentId: this.experimentId
      });
    },
    view() {
      return this.$store.getters["view/getView"]({
        viewId: this.viewId
      });
    }
  },
  methods: {
    refreshData() {
      this.errors = [];

      if (this.runId) {
        this.$store.dispatch("run/fetchRun", {runId: this.runId});
      }

      if (this.experimentId) {
        this.$store.dispatch("experiment/fetchExperiment", {experimentId: this.experimentId});
      }

      if (this.viewId) {
        this.$store.dispatch("view/fetchView", {viewId: this.viewId});
      }
    }
  },
  async mounted() {
    this.refreshData();
  }
};
</script>

<style scoped>

</style>
