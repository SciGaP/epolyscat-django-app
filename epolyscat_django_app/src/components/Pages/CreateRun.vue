<template>
  <div class="w-100 h-100 bg-light p-2">
    <div class="w-100 h-100 bg-white overflow-auto p-3">
      <b-breadcrumb>
        <template v-if="experimentIdFromQueryString">
          <b-breadcrumb-item to="/experiments">Experiments</b-breadcrumb-item>
          <b-breadcrumb-item v-if="experiment" :to="`/runs?experimentId=${experiment.experimentId}`">
            {{ experiment.name }}
          </b-breadcrumb-item>
        </template>
        <template v-else-if="viewId">
          <b-breadcrumb-item to="/Views">Views</b-breadcrumb-item>
          <b-breadcrumb-item v-if="view" :to="`/runs?experimentId=${view.viewId}`">
            {{ view.name }}
          </b-breadcrumb-item>
        </template>
        <b-breadcrumb-item v-else to="/runs">Runs</b-breadcrumb-item>

        <template v-if="cloneRunId">
          <b-breadcrumb-item :to="runLink" v-if="cloneRunOriginalName">
            {{ cloneRunOriginalName }}
          </b-breadcrumb-item>
          <b-breadcrumb-item :to="newRunLink">Clone</b-breadcrumb-item>
        </template>
        <b-breadcrumb-item v-else :to="newRunLink">New</b-breadcrumb-item>
      </b-breadcrumb>

      <div class="w-100">
        <div class="d-inline" style="font-weight: 400; font-size: 19px;">Create a new run</div>
      </div>


      <b-overlay :show="processing" class="w-100 d-flex flex-row">
        <div class="p-2 d-flex flex-row flex-fill">
          <div style="min-width: 150px;">
            <label for="root">Run root *</label>
          </div>
          <div>
            <b-overlay :show="cloneRunId && !cloneRun && !root" rounded="sm">
              <b-form-input
                  id="root"
                  size="sm"
                  v-model="root"
                  :state="inputState.root"
                  list="root-list"
                  autocomplete="off"
              />
            </b-overlay>
            <b-form-datalist id="root-list" :options="rootList"/>
            <b-form-invalid-feedback>
              The root of the run cannot be empty
            </b-form-invalid-feedback>
          </div>
        </div>
      </b-overlay>

      <RunViewableEditor :run-id="cloneRunId" :run="cloneRun" filename="inpc" :read-only="false"
                         v-on:change="onInpcContentChange"/>

      <b-overlay :show="processing || (cloneRunId && !cloneRun)" class="w-100">
        <div class="w-100 d-flex flex-row">
          <adpf-group-resource-profile-selector
              class="p-2 d-flex flex-row flex-fill"
              v-on:input="onGroupResourceProfileSelector"
              :value="groupResourceProfileId"
          />
          <adpf-experiment-compute-resource-selector
              v-if="trecxApplicationModuleId"
              class="p-2 d-flex flex-row flex-fill" v-on:input="onComputeResourceSelector"
              :application-module-id="trecxApplicationModuleId"
              :value="computeResourceId"
          />
        </div>
        <adpf-queue-settings-editor
            v-on:input="onQueueSettingEditor"
            :queue-name="queueName"
            :node-count="nodeCount"
            :total-cpu-count="coreCount"
            :wall-time-limit="wallTimeLimit"
            :total-physical-memory="totalPhysicalMemory"
        />
      </b-overlay>

      <div class="w-100 p-3 text-lg-left" style="font-size: 20px;">
        <b-form-invalid-feedback :state="inputState.root">
          <b-icon icon="info-circle-fill"/>&nbsp;
          Run root is required
        </b-form-invalid-feedback>
        <b-form-invalid-feedback :state="inputState.inpcContent">
          <b-icon icon="info-circle-fill"/>&nbsp;
          Run input is required
        </b-form-invalid-feedback>
        <b-form-invalid-feedback :state="inputState.computeResourceId">
          <b-icon icon="info-circle-fill"/>&nbsp;
          Compute resource is required
        </b-form-invalid-feedback>
        <b-form-invalid-feedback :state="inputState.groupResourceProfileId">
          <b-icon icon="info-circle-fill"/>&nbsp;
          Group resource profile is required
        </b-form-invalid-feedback>
        <b-form-invalid-feedback :state="inputState.queueName">
          <b-icon icon="info-circle-fill"/>&nbsp;
          A queue must be selected
        </b-form-invalid-feedback>
        <b-form-invalid-feedback :state="inputState.coreCount">
          <b-icon icon="info-circle-fill"/>&nbsp;
          Core count must be greater than zero
        </b-form-invalid-feedback>
        <b-form-invalid-feedback :state="inputState.nodeCount">
          <b-icon icon="info-circle-fill"/>&nbsp;
          Node count must greater than zero
        </b-form-invalid-feedback>
        <b-form-invalid-feedback :state="inputState.wallTimeLimit">
          <b-icon icon="info-circle-fill"/>&nbsp;
          Wall time limit must be specified
        </b-form-invalid-feedback>
        <b-form-invalid-feedback :state="inputState.totalPhysicalMemory">
          <b-icon icon="info-circle-fill"/>&nbsp;
          Total physical memory must be specified
        </b-form-invalid-feedback>
      </div>

      <div class="w-100 p-2 d-flex flex-row">
        <button-overlay :show="processing" class="w-100">
          <b-button variant="primary" class="w-100" v-on:click="onSave(false)">Save</b-button>
        </button-overlay>
        <button-overlay :show="processing" class="ml-2 w-100">
          <b-button variant="primary" class="w-100" v-on:click="onSave(true)">Submit</b-button>
        </button-overlay>
      </div>

    </div>
  </div>
</template>

<script>
import {ExperimentService, RunService} from "@/service/trecx-service";
import store from "@/store";
import ButtonOverlay from "@/components/overlay/button-overlay";
import RunViewableEditor from "@/components/block/RunViewableEditor";

export default {
  name: 'CreateRun',
  components: {RunViewableEditor, ButtonOverlay},
  store: store,
  data() {
    return {
      root: null,
      experimentId: null,
      groupResourceProfileId: null,
      computeResourceId: null,
      queueName: null,
      coreCount: null,
      nodeCount: null,
      wallTimeLimit: null,
      totalPhysicalMemory: null,

      inpcContent: null,
      inpcContentType: 'text',

      cloneRun: null,

      inputFieldsList: [
        "root",
        "inpcContent",
        "groupResourceProfileId",
        "computeResourceId",
        "queueName",
        "coreCount",
        "nodeCount",
        "wallTimeLimit",
        "totalPhysicalMemory",
      ],

      processing: false
    };
  },
  computed: {
    runLink() {
      let _link = "/runs/";
      if (this.cloneRunId) {
        _link += `${this.cloneRunId}?`;
      }

      if (this.experimentIdFromQueryString) {
        _link += `experimentId=${this.experimentIdFromQueryString}&`;
      }

      if (this.viewId) {
        _link += `viewId=${this.viewId}&`;
      }

      return _link;
    },
    inputState() {
      return {
        root: this.root === null ? null : this.isValid.root,
        inpcContent: this.inpcContent === null ? null : this.isValid.inpcContent,
        groupResourceProfileId: this.groupResourceProfileId === null ? null : this.isValid.groupResourceProfileId,
        computeResourceId: this.computeResourceId === null ? null : this.isValid.computeResourceId,
        queueName: this.queueName === null ? null : this.isValid.queueName,
        coreCount: this.coreCount === null ? null : this.isValid.coreCount,
        nodeCount: this.nodeCount === null ? null : this.isValid.nodeCount,
        wallTimeLimit: this.wallTimeLimit === null ? null : this.isValid.wallTimeLimit,
        totalPhysicalMemory: this.totalPhysicalMemory === null ? null : this.isValid.totalPhysicalMemory
      }
    },
    isValid() {
      return {
        root: !!this.root && this.root.length >= 1,
        inpcContent: this.inpcContentType !== "text" ? true : (!!this.inpcContent && this.inpcContent.length >= 1),
        groupResourceProfileId: !!this.groupResourceProfileId,
        computeResourceId: !!this.computeResourceId,
        queueName: !!this.queueName,
        coreCount: this.coreCount > 0,
        nodeCount: this.nodeCount > 0,
        wallTimeLimit: this.wallTimeLimit > 0,
        totalPhysicalMemory: this.totalPhysicalMemory >= 0
      }
    },
    isFormValid() {
      let _isFormValid = true;
      for (let i = 0; i < this.inputFieldsList.length; i++) {
        _isFormValid = _isFormValid && this.isValid[this.inputFieldsList[i]];
      }

      return _isFormValid;
    },
    newRunLink() {
      let newRunLink = "/create-run?";
      if (this.cloneRunId) {
        newRunLink += `cloneRunId=${this.cloneRunId}&`
      }
      if (this.experimentIdFromQueryString) {
        newRunLink += `experimentId=${this.experimentIdFromQueryString}&`
      }

      return newRunLink;
    },
    experimentIdFromQueryString() {
      return this.$route.query.experimentId;
    },
    viewId() {
      return this.$route.query.viewId;
    },
    cloneRunId() {
      return this.$route.query.cloneRunId;
    },
    cloneRunOriginalName() {
      const run = this.$store.getters["run/getRun"]({runId: this.cloneRunId});
      if (run) {
        return run.name
      } else {
        return null;
      }
    },
    experiments() {
      return this.$store.getters["experiment/getExperiments"]();
    },
    experimentsOptionsList() {
      if (this.experiments) {
        return this.experiments.map(({experimentId, name}) => {
          return {text: name, value: experimentId};
        });
      } else {
        return [];
      }
    },
    rootList() {
      if (this.experiments) {
        const roots = this.experiments.map(e => e.root);
        roots.sort();
        return roots;
      } else {
        return [];
      }
    },
    runs() {
      if (this.experimentId) {
        return this.$store.getters["run/getRuns"]({experimentId: this.experimentId});
      } else {
        return null;
      }
    },
    latestRunId() {
      if (this.cloneRunId) {
        return this.cloneRunId;
      } else if (this.runs && this.runs.length > 0) {
        return this.runs[this.runs.length - 1].runId;
      } else {
        return null;
      }
    },
    latestRun() {
      return this.$store.getters["run/getRun"]({runId: this.latestRunId});
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
    redirectLink(run) {
      let _link = "/runs?";

      if (run.experimentId) {
        _link += `experimentId=${run.experimentId}&`;
      }

      return _link;
    },
    makeFormVisited() {
      for (let i = 0; i < this.inputFieldsList.length; i++) {
        if (this[this.inputFieldsList[i]] === null) {
          this[this.inputFieldsList[i]] = "";
        }
      }
    },
    onGroupResourceProfileSelector(evt) {
      if (evt.detail && evt.detail.length > 0 && evt.detail[0]) {
        this.groupResourceProfileId = evt.detail[0];
      }
    },
    onComputeResourceSelector(evt) {
      if (evt.detail && evt.detail.length > 0 && evt.detail[0]) {
        this.computeResourceId = evt.detail[0];
      }
    },
    onQueueSettingEditor(evt) {
      if (evt.detail && evt.detail.length > 0 && evt.detail[0]) {
        this.queueName = evt.detail[0].queueName;
        this.coreCount = evt.detail[0].totalCPUCount;
        this.nodeCount = evt.detail[0].nodeCount;
        this.wallTimeLimit = evt.detail[0].wallTimeLimit;
        this.totalPhysicalMemory = evt.detail[0].totalPhysicalMemory;
      }
    },
    async onSave(submit = false) {
      this.makeFormVisited();
      if (this.isFormValid) {
        this.processing = true;

        const createRunPayload = {
          root: this.root,
          experimentId: this.experimentId,
          groupResourceProfileId: this.groupResourceProfileId,
          computeResourceId: this.computeResourceId,
          queueName: this.queueName,
          coreCount: this.coreCount,
          nodeCount: this.nodeCount,
          wallTimeLimit: this.wallTimeLimit,
          totalPhysicalMemory: this.totalPhysicalMemory
        };

        if (this.inpcContentType === "text") {
          createRunPayload.directedit = this.inpcContent;
        } else if (this.inpcContentType === "table") {
          createRunPayload.inputTable = this.inpcContent;
        }

        const run = await RunService.createRun(createRunPayload, submit);
        if (!submit) {
          // Refresh the views so that "Unsubmitted" list is refreshed.
          this.$store.dispatch("view/fetchViews");
        }

        // If creating run resulted in the creation of a new root/experiment,
        // reload the list of experiments
        if (this.experimentId !== run.experimentId) {
          this.$store.dispatch("experiment/fetchExperiments");
        }
        this.$router.history.push(this.redirectLink(run));

        this.processing = false;
      }
    },
    refreshData() {
      ExperimentService.fetchAllExperiments();

      if (this.experimentId) {
        this.$store.dispatch("run/fetchRuns", {experimentId: this.experimentId});
      }

      if (this.viewId) {
        this.$store.dispatch("view/fetchView", {viewId: this.viewId});
      }

      this.refreshRun();
    },
    async refreshRun() {
      if (this.latestRunId) {
        this.$store.dispatch("run/fetchRun", {runId: this.latestRunId});

        this.cloneRun = await RunService.cloneRun({runId: this.latestRunId});
        this.experimentId = this.cloneRun.experimentId;
      }
    },
    onInpcContentChange(inpcContentType, inpcContent) {
      this.inpcContentType = inpcContentType;
      this.inpcContent = inpcContent;
    }
  },
  watch: {
    async experiment() {
      if (this.experiment) this.root = this.experiment.root;
    },
    latestRunId() {
      this.refreshRun();
    },
    experimentIdFromQueryString() {
      this.experimentId = this.experimentIdFromQueryString;
    },
    latestRun() {
      // Copy compute resource details if they exist
      if (this.latestRun && this.latestRun.groupResourceProfileId && this.latestRun.computeResourceId && this.latestRun.queueName) {
        this.groupResourceProfileId = this.latestRun.groupResourceProfileId;
        this.computeResourceId = this.latestRun.computeResourceId;
        this.queueName = this.latestRun.queueName;
        this.coreCount = this.latestRun.coreCount;
        this.nodeCount = this.latestRun.nodeCount;
        this.wallTimeLimit = this.latestRun.wallTimeLimit;
        this.totalPhysicalMemory = this.latestRun.totalPhysicalMemory;
      }
    }
  },
  async mounted() {
    this.experimentId = this.experimentIdFromQueryString;

    this.refreshData();

    this.$store.dispatch("settings/fetchSettings");
  }
};
</script>

<style scoped>

.visible {
  visibility: unset;
}

.invisible {
  visibility: hidden;
  position: fixed;
  top: -10000px;
}

/*.run-pane {*/
/*  width: 100%;*/
/*  height: 300px;*/
/*  padding-top: 30px;*/
/*}*/

/*.create-run {*/
/*  background-color: #FFFFFF;*/
/*  border-radius: 10px;*/
/*  height: 100%;*/
/*  overflow-y: scroll;*/
/*  padding: 48px 80px 50px 50px;*/
/*}*/

/*input[type=text] {*/
/*  border-radius: 5px;*/
/*  width: 100%;*/
/*  height: 47px;*/
/*  border: 0.5px solid #838383;*/
/*  background: #FFFFFF;*/
/*  text-align: inherit;*/
/*}*/

/*.filter-badge {*/
/*  background-color: #F6F6F6;*/
/*  border: 1px solid #F6F6F6;*/
/*  border-radius: 30px;*/
/*  text-align: center;*/
/*  padding: 10px 20px;*/
/*  color: #000000;*/
/*}*/

/*.run-submit {*/
/*  width: 100%;*/
/*  text-align: center;*/
/*  padding-top: 45px !important;*/
/*}*/

/*.run-submit-btn {*/
/*  border: 1px solid #226597;*/
/*  box-sizing: border-box;*/
/*  border-radius: 10px;*/
/*  width: 100%;*/
/*  text-align: center;*/
/*  height: 55px !important;*/
/*  vertical-align: center;*/
/*}*/

/*.active {*/
/*  background-color: #226597;*/
/*  color: #FFFFFF;*/
/*}*/

/*.run-dropdown::after {*/
/*  margin-left: 95%;*/
/*}*/

/*.grid-container--2 {*/
/*  width: 100%;*/
/*  column-gap: 58px;*/
/*  row-gap: 30px;*/
/*  grid-template-columns: 1fr 1fr !important;*/
/*  display: grid;*/
/*  padding-top: 53px !important;*/
/*}*/

/*.run-dropdown {*/
/*  width: 100%;*/
/*  background: #FFFFFF;*/
/*  border: 0.5px solid #838383;*/
/*  box-sizing: border-box;*/
/*  border-radius: 5px;*/
/*  height: 47px;*/
/*}*/

/*.check-settings {*/
/*  background: #226597;*/
/*  box-shadow: 0px 4px 20px rgba(208, 208, 208, 0.2);*/
/*  border-radius: 5px;*/
/*  width: 37px;*/
/*  height: 35px;*/
/*  padding: 4px 10px;*/
/*}*/

/*.settings {*/
/*  background-color: #F6F6F6;*/
/*  border-radius: 5px;*/
/*  padding: 13px 30px 26px 20px;*/
/*  box-sizing: border-box;*/
/*}*/

/*.settings-grid {*/
/*  padding-top: 20px;*/
/*  display: grid;*/
/*  column-gap: 20px !important;*/
/*  grid-template-columns: 1fr 1fr 1fr 1fr !important;*/
/*}*/

/*.grid-box-settings {*/
/*  width: 100%;*/
/*  text-align: center;*/
/*}*/

/*input:disabled {*/
/*  border: none;*/
/*  background-color: #F6F6F6;*/

/*}*/

/*.custom-resizer {*/
/*  width: 100%;*/
/*  height: 100%;*/
/*}*/

/*.custom-resizer > .multipane-resizer {*/
/*  margin: 0;*/
/*  left: 0;*/
/*  position: relative;*/
/*}*/

/*.custom-resizer > .multipane-resizer:hover::before {*/
/*  border-color: #237ABB;*/
/*}*/

/*.custom-resizer > .multipane-resizer::before {*/
/*  display: block;*/
/*  content: "";*/
/*  width: 9px;*/
/*  height: 67px;*/
/*  position: absolute;*/
/*  top: 50%;*/
/*  left: 50%;*/
/*  margin-top: -30px;*/
/*  margin-left: -10px;*/
/*  border-left: 10px solid #237ABB;*/
/*  border-right: 10px solid #237ABB;*/
/*  border-radius: 20px;*/
/*  z-index: 1;*/
/*}*/

/*.multipane > div {*/
/*  z-index: auto;*/
/*}*/
</style>
