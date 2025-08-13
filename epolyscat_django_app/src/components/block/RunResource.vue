<template>
  <div class="w-100">
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
        :disabled="true"
        v-on:input="onQueueSettingEditor"
        :queue-name="queueName"
        :node-count="nodeCount"
        :total-cpu-count="coreCount"
        :wall-time-limit="wallTimeLimit"
        :total-physical-memory="totalPhysicalMemory"
    />
  </div>
</template>

<script>
import store from "@/store";

export default {
  name: 'RunResource',
  props: {
    run: {},
    disabled: {
      default: false
    }
  },
  store: store,
  data() {
    return {
      groupResourceProfileId: null,
      computeResourceId: null,
      queueName: null,
      coreCount: null,
      nodeCount: null,
      wallTimeLimit: null,
      totalPhysicalMemory: null
    };
  },
  computed: {
    isResourceSelectionAllowed() {
      return !this.disabled && !!this.run && ["Unsubmitted", "FAILED"].indexOf(this.run.status) >= 0;
    },
    trecxApplicationModuleId() {
      return this.$store.getters["settings/trecxApplicationModuleId"];
    }
  },
  methods: {
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
    updateRun() {
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
      this.updateRun();
    }
  },
  async beforeMount() {
    this.updateRun();
    this.$store.dispatch("settings/fetchSettings");
  }
};
</script>

<style scoped>
</style>
