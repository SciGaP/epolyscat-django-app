<template>
  <div class="w-100 h-100 bg-light p-2">
    <multipane class="w-100 h-100" layout="vertical">
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
          </b-breadcrumb>

          <Errors :errors="errors"/>

        </div>
        <div class="w-100 overflow-auto">
          <div class="w-100 d-flex mb-3">
            <div class="d-inline pr-3 flex-fill overflow-auto" style="font-weight: 400; font-size: 19px;">
              <span v-if="run">{{ run.name }}</span>
            </div>
            <RunActions :run="run" :experiment="experiment" :view="view" variant="outline-primary"
                        :showButtonText="true" v-on:delete="onDelete"/>
          </div>

          <div class="w-100 d-flex border" style="border-radius: 10px;">
            <div class="flex-fill p-2">
              <div class="text-center text-primary pb-2" style="font-weight: 600;">Viewables</div>
              <b-overlay :show="!viewables">
                <div class="text-center" style="min-height: 50px;">
                  <b-button pill :variant="viewable.filename === filename ? 'primary' : 'outline-primary'" class="mb-2"
                            v-for="(viewable, viewableIndex) in viewables"
                            :key="viewableIndex" v-on:click="filename = viewable.filename">
                    {{ viewable.filename }}
                  </b-button>
                </div>
              </b-overlay>
            </div>
            <div class="flex-fill p-2 border-left border-right">
              <div class="text-center text-primary pb-2" style="font-weight: 600;">Input Files</div>
              <b-overlay :show="!inputFiles">
                <div class="text-center" style="min-height: 50px;">
                  <b-button pill variant="outline-primary" v-for="(inputFile, inputFileIndex) in inputFiles"
                            :key="inputFileIndex" class="mb-2">
                    {{ inputFile.filename }}
                  </b-button>
                </div>
              </b-overlay>
            </div>
            <div class="flex-fill p-2">
              <div class="text-center text-primary pb-2" style="font-weight: 600;">Plotables</div>
              <b-overlay :show="!plotables">
                <div class="text-center" style="min-height: 50px;">
                  <b-button pill :variant="selectedPlotable === plotable ? 'primary' : 'outline-primary'"
                            v-for="(plotable, plotableIndex) in plotables" class="mb-2"
                            :key="plotableIndex" v-on:click="togglePlotableSelection(plotable)">
                    {{ plotable }}
                  </b-button>
                </div>
              </b-overlay>
            </div>
          </div>

          <RunViewableEditor :run-id="runId" :run="run" :filename="filename" :filenames="filenames" :read-only="true"
                             :allow-refresh="true"/>

          <RunResource :run="run" :disabled="true"/>
        </div>
      </div>
      <multipane-resizer v-if="selectedPlotable"/>
      <div class="h-100 m-1" v-if="selectedPlotable">
        <PostProcessing :selectedRunIds="[runId]" :plotable="selectedPlotable"
                        v-on:plotable-change="setPlotableSelected"/>
      </div>
    </multipane>
  </div>
</template>

<script>
import store from "@/store";
import RunActions from "@/components/block/RunActions";
import {PlotService} from "@/service/trecx-service";
import Errors from "@/components/block/Errors";
import {Multipane, MultipaneResizer} from "vue-multipane";
import PostProcessing from "@/components/block/PostProcessing";
import RunViewableEditor from "@/components/block/RunViewableEditor";
import RunResource from "@/components/block/RunResource";

export default {
  name: 'CreateRun',
  components: {RunResource, RunViewableEditor, PostProcessing, Errors, RunActions, Multipane, MultipaneResizer},
  store: store,
  data() {
    return {
      plotables: null,
      viewables: null,
      filename: "inpc",
      inputFiles: null,
      processingPlotables: null,
      processingViewables: null,
      processingInputFiles: null,
      processing: false,

      selectedPlotable: null,

      errors: [],

      runRefreshTimeout: null
    };
  },
  computed: {
    filenames() {
      if (this.viewables) {
        return this.viewables.map(({filename}) => filename);
      } else {
        return null;
      }
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
    },
    trecxApplicationModuleId() {
      return this.$store.getters["settings/trecxApplicationModuleId"];
    },
  },
  methods: {
    viewableLink({filename}) {
      let _link = `/runs/${this.runId}/viewables/${filename}?`;

      if (this.experimentId) {
        _link += `experimentId=${this.experimentId}&`;
      }

      if (this.viewId) {
        _link += `viewId=${this.viewId}&`;
      }

      return _link;
    },
    onDelete() {
      let redirectUrl = "/runs?";
      if (this.experimentId) {
        redirectUrl += `experimentId=${this.experimentId}&`;
      }

      if (this.viewId) {
        redirectUrl += `viewId=${this.viewId}&`;
      }

      this.$router.history.push(redirectUrl);
    },
    refreshData() {
      this.errors = [];

      if (this.runId) {
        this.$store.dispatch("run/fetchRun", {runId: this.runId});

        PlotService.getPlotables({runIds: [this.runId]}).then(plotables => {
          this.plotables = plotables;
          this.processingPlotables = false;
        }).catch(e => {
          this.processingPlotables = false;
          this.errors.push({
            variant: "danger",
            title: "Network Error",
            description: "Unknown error when fetching plotables.",
            source: e
          });
        });

        PlotService.getViewables({runId: this.runId}).then(viewables => {
          this.viewables = viewables;
          this.processingViewables = false;
        }).catch(e => {
          this.processingViewables = false;
          this.errors.push({
            variant: "danger",
            title: "Network Error",
            description: "Unknown error when fetching viewables.",
            source: e
          });
        });

        PlotService.getInputFiles({runId: this.runId}).then(inputFiles => {
          this.inputFiles = inputFiles;
          this.processingInputFiles = false;
        }).catch(e => {
          this.processingInputFiles = false;
          this.errors.push({
            variant: "danger",
            title: "Network Error",
            description: "Unknown error when fetching input files.",
            source: e
          });
        });
      }

      if (this.experimentId) {
        this.$store.dispatch("experiment/fetchExperiment", {experimentId: this.experimentId});
      }

      if (this.viewId) {
        this.$store.dispatch("view/fetchView", {viewId: this.viewId});
      }
    },
    setPlotableSelected(plotable) {
      this.selectedPlotable = plotable;
    },
    togglePlotableSelection(plotable) {
      if (this.selectedPlotable !== plotable) {
        this.setPlotableSelected(plotable);
      } else {
        // Toggle the selection.
        this.selectedPlotable = null;
      }
    }
  },
  async mounted() {
    this.errors = [];

    this.plotables = null;
    this.processingPlotables = true;

    this.viewables = null;
    this.processingViewables = true;

    this.inputFiles = null;
    this.processingInputFiles = true;

    this.refreshData();

    this.runRefreshTimeout = setInterval(() => {
      this.refreshData();
    }, 15000);
  },
  beforeDestroy() {
    clearInterval(this.runRefreshTimeout);
  }
};
</script>

<style scoped>
.multipane .multipane-resizer {
  margin: 0;
  left: 0;
  position: relative;
}

.multipane .multipane-resizer:before {
  display: block;
  content: "";
  width: 90%;
  height: 40px;
  position: absolute;
  top: 50%;
  left: 5%;
  margin-top: -20px;
  margin-left: -10%;
  background-color: #237abb;
  border-radius: 10px;
}

.multipane .multipane-resizer:hover,
.multipane .multipane-resizer:before {
  border-color: #999;
}
</style>
