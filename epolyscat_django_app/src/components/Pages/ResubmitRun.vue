<template>
  <div class="w-100 h-100 bg-light p-2">
    <div class="w-100 h-100 bg-white overflow-auto p-3">
      <b-breadcrumb>
        <template v-if="experimentId">
          <b-breadcrumb-item to="/experiments">Experiments</b-breadcrumb-item>
          <b-breadcrumb-item v-if="experiment" :to="`/runs?experimentId=${experiment.experimentId}`">
            {{ experiment.name }}
          </b-breadcrumb-item>
        </template>
        <b-breadcrumb-item v-else to="/runs">Runs</b-breadcrumb-item>

        <b-breadcrumb-item :to="`/runs/${runId}`" v-if="run">{{ run.name }}</b-breadcrumb-item>
        <b-breadcrumb-item :to="`/resubmit-run?runId=${runId}`" v-if="run">
          <template v-if="run.canSubmit">Submit</template>
          <template v-else>Resubmit</template>
        </b-breadcrumb-item>
      </b-breadcrumb>

      <div class="w-100">
        <div class="d-inline" style="font-weight: 400; font-size: 19px;" v-if="run">
          <template v-if="run.canSubmit">Submit run</template>
          <template v-else>Resubmit run</template>
        </div>
      </div>

      <div class="w-100 d-flex flex-row">
        <div class="p-2 d-flex flex-row flex-fill">
          <div class="p-2 d-flex flex-row flex-fill">
            <div class="pr-3"><label for="name">Name</label></div>
            <div id="name" v-if="run">
              <b-form-input
                  id="root"
                  size="sm"
                  v-model="run.name"
                  :disabled="true"
              />
            </div>
          </div>
          <div class="p-2 d-flex flex-row flex-fill">
            <div class="pr-3"><label for="status">Status</label></div>
            <div id="status" v-if="run">
              <b-form-input
                  id="root"
                  size="sm"
                  v-model="run.status"
                  :disabled="true"
              />
            </div>
          </div>
        </div>
      </div>

      <RunViewableEditor :run-id="runId" :run="run" filename="inpc"/>

      <b-overlay :show="processing || !run" class="w-100">
        <div class="w-100 d-flex flex-row">
          <adpf-group-resource-profile-selector
              class="p-2 d-flex flex-row flex-fill" :disabled="!isResourceSelectionAllowed"
              v-on:input="onGroupResourceProfileSelector"
              :value="groupResourceProfileId"
          />
          <adpf-experiment-compute-resource-selector
              v-if="trecxApplicationModuleId" :disabled="!isResourceSelectionAllowed"
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

      <div class="w-100 p-2">
        <button-overlay :show="processing" class="w-100">
          <b-button variant="primary" class="w-100" v-on:click="onSubmit">Submit</b-button>
        </button-overlay>
      </div>

    </div>
  </div>
</template>

<script>
import {RunService} from "@/service/trecx-service";
import store from "@/store";
import ButtonOverlay from "@/components/overlay/button-overlay";
import RunViewableEditor from "@/components/block/RunViewableEditor";

export default {
  name: 'CreateRun',
  components: {RunViewableEditor, ButtonOverlay},
  store: store,
  data() {
    return {
      groupResourceProfileId: null,
      computeResourceId: null,
      queueName: null,
      coreCount: null,
      nodeCount: null,
      wallTimeLimit: null,
      totalPhysicalMemory: null,

      inputFieldsList: [
        "groupResourceProfileId",
        "computeResourceId",
        "queueName",
        "coreCount",
        "nodeCount",
        "wallTimeLimit",
        "totalPhysicalMemory",
      ],

      processing: false,
    };
  },
  computed: {
    isResourceSelectionAllowed() {
      return !!this.run && ["Unsubmitted", "FAILED"].indexOf(this.run.status) >= 0;
    },
    inputState() {
      return {
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
    runId() {
      return this.$route.query.runId;
    },
    experimentId() {
      return this.$route.query.experimentId;
    },
    run() {
      if (this.runId) {
        return this.$store.getters["run/getRun"]({runId: this.runId});
      } else {
        return null;
      }
    },
    experiment() {
      if (this.experimentId) {
        return this.$store.getters["experiment/getExperiment"]({experimentId: this.experimentId});
      } else {
        return null;
      }
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
    async onSubmit() {
      this.makeFormVisited();
      if (this.isFormValid) {
        this.processing = true;

        const resubmitRunPayload = {
          runId: this.runId,
          groupResourceProfileId: this.groupResourceProfileId,
          computeResourceId: this.computeResourceId,
          queueName: this.queueName,
          coreCount: this.coreCount,
          nodeCount: this.nodeCount,
          wallTimeLimit: this.wallTimeLimit,
          totalPhysicalMemory: this.totalPhysicalMemory
        };

        if (this.run.canSubmit) {
          await RunService.submitRun(resubmitRunPayload);
        } else {
          await RunService.resubmitRun(resubmitRunPayload);
        }

        this.$router.history.push(this.redirectLink(this.run));

        this.processing = false;
      }
    },
    refreshData() {
      this.onRunValueChange();

      if (this.runId) {
        this.$store.dispatch("run/fetchRun", {runId: this.runId});
      }

      if (this.experimentId) {
        this.$store.dispatch("experiment/fetchExperiment", {experimentId: this.experimentId});
      }
    },
    onRunValueChange() {
      if (this.run) {
        this.groupResourceProfileId = this.run.groupResourceProfileId;
        this.computeResourceId = this.run.computeResourceId;
        this.nodeCount = this.run.nodeCount;
        this.queueName = this.run.queueName;
        this.coreCount = this.run.coreCount;
        this.wallTimeLimit = this.run.wallTimeLimit;
        this.totalPhysicalMemory = this.run.totalPhysicalMemory;
      }
    }
  },
  watch: {
    run() {
      this.onRunValueChange();
    }
  },
  async mounted() {
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
</style>
