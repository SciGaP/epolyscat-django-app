<template>
  <div>
    <div class="w-100 d-flex flex-row">
      <adpf-group-resource-profile-selector class="p-2 d-flex flex-row flex-fill" :disabled="viewing"
        v-on:input="onGroupResourceProfileSelector" :value="groupResourceProfileId" />
      <adpf-experiment-compute-resource-selector v-if="epolyscatApplicationModuleId" :disabled="viewing"
        class="p-2 d-flex flex-row flex-fill" v-on:input="onComputeResourceSelector"
        :application-module-id="epolyscatApplicationModuleId" :value="computeResourceId" />
    </div>
    <adpf-queue-settings-editor :disabled="viewing" v-on:input="onQueueSettingEditor" :queue-name="queueName"
      :node-count="nodeCount" :total-cpu-count="coreCount" :wall-time-limit="wallTimeLimit"
      :total-physical-memory="totalPhysicalMemory" :style="{'pointer-events': (viewing) ? 'none' : 'initial'}"/>
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
        },
        viewing: Boolean
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
            return !this.disabled && !!this.run && ["UNSUBMITTED", "FAILED"].indexOf(this.run.status) >= 0;
        },
        epolyscatApplicationModuleId() {
            return this.$store.getters["settings/epolyscatApplicationModuleId"];
        }
    },
    methods: {
        onGroupResourceProfileSelector(evt) {
            if (evt.detail && evt.detail.length > 0 && evt.detail[0] && !this.viewing) {
                this.groupResourceProfileId = evt.detail[0];
            }
            // console.log("FROM: onGroupResourceProfileSelector")

            this.emitResources();
        },
        onComputeResourceSelector(evt) {
            if (evt.detail && evt.detail.length > 0 && evt.detail[0] && !this.viewing) {
                this.computeResourceId = evt.detail[0];
            }
            // console.log("FROM: onComputeResourceSelector")

            this.emitResources();
        },
        onQueueSettingEditor(evt) {
            // console.log(this.numTimesUpdated, " -- FROM: onQueueSettingEditor", evt);

            if (evt.detail && evt.detail.length > 0 && evt.detail[0] && !this.viewing) {
                this.queueName = evt.detail[0].queueName;
                this.coreCount = evt.detail[0].totalCPUCount;
                this.nodeCount = evt.detail[0].nodeCount;
                this.wallTimeLimit = evt.detail[0].wallTimeLimit;
                this.totalPhysicalMemory = evt.detail[0].totalPhysicalMemory;
            }
                
            this.emitResources();
        },
        emitResources() {
            let updatedResources = {
                groupResourceProfileId: this.groupResourceProfileId,
                computeResourceId: this.computeResourceId,
                nodeCount: parseInt(this.nodeCount),
                queueName: this.queueName,
                coreCount: parseInt(this.coreCount),
                wallTimeLimit: parseInt(this.wallTimeLimit),
                totalPhysicalMemory: parseInt(this.totalPhysicalMemory),
            };

            // console.log("RESOURCES UPDATING!!", updatedResources, "run previously", {...this.run}, "this.viewing:", this.viewing);
            if (!this.viewing)
                this.$emit("updateResources", updatedResources);
        },
        updateFromRun() {
            if (this.run) {
                this.groupResourceProfileId = this.run.groupResourceProfileId;
                this.computeResourceId = this.run.computeResourceId;
                this.nodeCount = this.run.nodeCount;
                this.queueName = this.run.queueName;
                this.coreCount = this.run.coreCount;
                this.wallTimeLimit = this.run.wallTimeLimit;
                this.totalPhysicalMemory = this.run.totalPhysicalMemory;
            }
            // console.log("updateFromRun (after):", {...this}, "run: ", {...this.run});
        }
    },
    watch: {
        run() {
            this.updateFromRun();
        },
        queueName(newer, older) {
            // console.log("queueName, new:", newer, "old:", older)
        }
    },
    async beforeMount() {
        this.updateFromRun();
        await this.$store.dispatch("settings/fetchSettings");
    }
};
</script>
