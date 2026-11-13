<template>
  <section class="resource-settings-grid">
    <div class="resource-fields">
      <div class="resource-row">
        <label for="allocation">Allocation</label>
        <div class="resource-search-dropdown" :class="{ open: showAllocationMenu }">
          <b-form-input
              id="allocation"
              v-model="allocationSearch"
              class="resource-search-input"
              :state="resourceInputState.groupResourceProfileId"
              placeholder="Select"
              autocomplete="off"
              :disabled="disabled"
              v-on:focus="openResourceMenu('allocation')"
              v-on:input="onAllocationSearchInput"
          />
          <button
              type="button"
              class="resource-search-toggle"
              :disabled="disabled"
              v-on:click="toggleResourceMenu('allocation')"
          >
            <b-icon icon="chevron-down" />
          </button>
          <div v-if="showAllocationMenu" class="resource-search-menu">
            <b-form-input
                v-model="allocationSearch"
                class="resource-menu-search"
                placeholder="Search allocation"
                v-on:click.stop
                v-on:input="onAllocationSearchInput"
            />
            <button
                v-for="option in filteredAllocationOptions"
                :key="option.value"
                type="button"
                class="resource-option"
                :class="{ active: option.value === groupResourceProfileId }"
                v-on:click="selectAllocationOption(option)"
            >
              <span class="resource-option-label">{{ option.label }}</span>
            </button>
            <div v-if="!filteredAllocationOptions.length" class="resource-empty-state">
              Select an available option
            </div>
          </div>
          <b-form-invalid-feedback :state="resourceInputState.groupResourceProfileId">
            Group resource profile is required
          </b-form-invalid-feedback>
        </div>
      </div>

      <div v-if="showSaveTarget" class="resource-row">
        <label for="save-to">Save to</label>
        <div class="save-to-dropdown" :class="{ open: showSaveTargetMenu }">
          <button
              id="save-to"
              type="button"
              class="save-to-trigger"
              :disabled="disabled"
              v-on:click="showSaveTargetMenu = !showSaveTargetMenu"
          >
            <span>{{ selectedSaveTargetLabel }}</span>
            <b-icon icon="chevron-down" />
          </button>
          <div v-if="showSaveTargetMenu" class="save-to-menu">
            <b-form-input
                v-model="saveTargetSearch"
                class="save-to-search"
                placeholder="Search here"
                v-on:click.stop
            />
            <button
                v-for="option in filteredSaveTargetOptions"
                :key="option"
                type="button"
                class="save-to-option"
                v-on:click="selectSaveTarget(option)"
            >
              {{ option }}
            </button>
          </div>
        </div>
      </div>
    </div>

    <div class="resource-fields">
      <div class="resource-row">
        <label for="compute-resource">Compute Resource</label>
        <div class="resource-search-dropdown" :class="{ open: showComputeResourceMenu }">
          <b-form-input
              id="compute-resource"
              v-model="computeResourceSearch"
              class="resource-search-input"
              :state="resourceInputState.computeResourceId"
              placeholder="Select"
              autocomplete="off"
              :disabled="disabled"
              v-on:focus="openResourceMenu('compute')"
              v-on:input="onComputeResourceSearchInput"
          />
          <button
              type="button"
              class="resource-search-toggle"
              :disabled="disabled"
              v-on:click="toggleResourceMenu('compute')"
          >
            <b-icon icon="chevron-down" />
          </button>
          <div v-if="showComputeResourceMenu" class="resource-search-menu">
            <b-form-input
                v-model="computeResourceSearch"
                class="resource-menu-search"
                placeholder="Search compute resource"
                v-on:click.stop
                v-on:input="onComputeResourceSearchInput"
            />
            <button
                v-for="option in filteredComputeResourceOptions"
                :key="option.value"
                type="button"
                class="resource-option"
                :class="{ active: option.value === computeResourceId }"
                v-on:click="selectComputeResourceOption(option)"
            >
              <span class="resource-option-label">{{ option.label }}</span>
            </button>
            <div v-if="!filteredComputeResourceOptions.length" class="resource-empty-state">
              {{ computeResourceEmptyMessage }}
            </div>
          </div>
          <b-form-invalid-feedback :state="resourceInputState.computeResourceId">
            {{ computeResourceValidationMessage }}
          </b-form-invalid-feedback>
          <div
              class="resource-readiness"
              :class="`status-${resourceReadiness.status}`"
              role="status"
              aria-live="polite"
          >
            <b-icon :icon="resourceReadinessIcon" aria-hidden="true" />
            <span>{{ resourceReadiness.message }}</span>
          </div>
        </div>
      </div>
    </div>

    <div class="settings-label">Settings</div>
    <div class="queue-settings-card">
      <div class="queue-card-header">
        <div class="queue-name">{{ queueDisplayName }}</div>
        <div class="queue-selector" :class="{ open: showQueueMenu }">
          <button
              type="button"
              class="queue-check"
              :class="{ active: showQueueMenu }"
              :disabled="disabled || queueOptions.length === 0"
              title="Select queue"
              v-on:click="toggleQueueMenu"
          >
            <b-icon icon="chevron-down" />
          </button>
          <div v-if="showQueueMenu" class="queue-menu">
            <b-form-input
                v-model="queueSearch"
                class="queue-search"
                placeholder="Search queue"
                v-on:click.stop
            />
            <button
                v-for="option in filteredQueueOptions"
                :key="option.value"
                type="button"
                class="queue-option"
                :class="{ active: option.value === queueName }"
                v-on:click="selectQueueOption(option)"
            >
              {{ option.label }}
            </button>
            <div v-if="!filteredQueueOptions.length" class="queue-empty-state">
              Select an available option
            </div>
          </div>
        </div>
      </div>
      <div class="queue-fields">
        <label>
          <b-form-input
              v-model.number="coreCount"
              type="number"
              class="queue-input"
              :state="resourceInputState.coreCount"
              :disabled="disabled"
              v-on:input="onCoreCountInput"
          />
          <span>Core count</span>
        </label>
        <button
            type="button"
            class="queue-lock-button"
            :class="{ active: queueResourceLockEnabled }"
            :disabled="disabled || !selectedQueueDefault || selectedQueueDefault.cpuPerNode <= 0"
            title="Link node and core counts"
            v-on:click="toggleQueueResourceLock"
        >
          <b-icon :icon="queueResourceLockEnabled ? 'lock' : 'unlock'" />
        </button>
        <label>
          <b-form-input
              v-model.number="nodeCount"
              type="number"
              class="queue-input"
              :state="resourceInputState.nodeCount"
              :disabled="disabled"
              v-on:input="onNodeCountInput"
          />
          <span>Node count</span>
        </label>
        <label>
          <b-form-input
              v-model.number="wallTimeLimit"
              type="number"
              class="queue-input"
              :state="resourceInputState.wallTimeLimit"
              :disabled="disabled"
          />
          <span>Wall time limit</span>
        </label>
        <label>
          <b-form-input
              v-model="totalPhysicalMemory"
              class="queue-input memory-input"
              :state="resourceInputState.totalPhysicalMemory"
              :disabled="disabled"
          />
          <span>Total physical memory</span>
        </label>
      </div>
    </div>

    <div class="resource-validation-copy" aria-hidden="true">
      <b-form-invalid-feedback :state="resourceInputState.queueName">
        A queue must be selected
      </b-form-invalid-feedback>
      <b-form-invalid-feedback :state="resourceInputState.coreCount">
        Core count must be greater than zero
      </b-form-invalid-feedback>
      <b-form-invalid-feedback :state="resourceInputState.nodeCount">
        Node count must greater than zero
      </b-form-invalid-feedback>
      <b-form-invalid-feedback :state="resourceInputState.wallTimeLimit">
        Wall time limit must be specified
      </b-form-invalid-feedback>
      <b-form-invalid-feedback :state="resourceInputState.totalPhysicalMemory">
        Total physical memory must be specified
      </b-form-invalid-feedback>
    </div>

    <div class="native-resource-selectors" aria-hidden="true">
      <adpf-queue-settings-editor
          :queue-name="queueName"
          :node-count="nodeCount"
          :total-cpu-count="coreCount"
          :wall-time-limit="wallTimeLimit"
          :total-physical-memory="totalPhysicalMemory"
      />
    </div>
  </section>
</template>

<script>
import store from "@/store";
import {
  buildComputeResourceReadiness,
  filterEligibleApplicationDeployments,
} from "@/utils/compute-resource-readiness";

const {services} = AiravataAPI;

export default {
  name: "RunResourceSettings",
  store,
  props: {
    run: {
      type: Object,
      required: true,
    },
    inputState: {
      type: Object,
      default: () => ({}),
    },
    applicationModuleId: {
      type: String,
      default: null,
    },
    applicationLabel: {
      type: String,
      default: "Application",
    },
    deploymentExecutionKind: {
      type: String,
      default: "module",
      validator: value => ["module", "utility"].includes(value),
    },
    disabled: {
      type: Boolean,
      default: false,
    },
    showSaveTarget: {
      type: Boolean,
      default: false,
    },
    saveTarget: {
      type: String,
      default: null,
    },
    saveTargetOptions: {
      type: Array,
      default: () => [],
    },
  },
  data() {
    return {
      allocationSearch: "",
      computeResourceSearch: "",
      showAllocationMenu: false,
      showComputeResourceMenu: false,
      groupResourceProfiles: [],
      applicationDeployments: [],
      allocationOptions: [],
      computeResourceOptions: [],
      queueOptions: [],
      queueSearch: "",
      showQueueMenu: false,
      queueOptionsLoading: false,
      queueResourceLockEnabled: true,
      resourceOptionsLoading: false,
      resourceOptionsLoaded: false,
      resourceOptionsError: null,
      computeResourcesLoading: false,
      computeResourceLoadError: null,
      computeResourceRequestId: 0,
      loadedApplicationModuleId: null,
      loadedGroupResourceProfileId: null,
      loadedDeploymentExecutionKind: null,
      saveTargetSearch: "",
      showSaveTargetMenu: false,
    };
  },
  computed: {
    groupResourceProfileId: {
      get() {
        return this.run ? this.run.groupResourceProfileId : null;
      },
      set(value) {
        this.updateRunResources({ groupResourceProfileId: value });
      }
    },
    computeResourceId: {
      get() {
        return this.run ? this.run.computeResourceId : null;
      },
      set(value) {
        this.updateRunResources({ computeResourceId: value });
      }
    },
    queueName: {
      get() {
        return this.run ? this.run.queueName : null;
      },
      set(value) {
        this.updateRunResources({ queueName: value });
      }
    },
    coreCount: {
      get() {
        return this.run ? this.run.coreCount : null;
      },
      set(value) {
        this.updateRunResources({ coreCount: value });
      }
    },
    nodeCount: {
      get() {
        return this.run ? this.run.nodeCount : null;
      },
      set(value) {
        this.updateRunResources({ nodeCount: value });
      }
    },
    wallTimeLimit: {
      get() {
        return this.run ? this.run.wallTimeLimit : null;
      },
      set(value) {
        this.updateRunResources({ wallTimeLimit: value });
      }
    },
    totalPhysicalMemory: {
      get() {
        return this.run ? this.run.totalPhysicalMemory : null;
      },
      set(value) {
        this.updateRunResources({ totalPhysicalMemory: value });
      }
    },
    resourceInputState() {
      return {
        groupResourceProfileId: null,
        computeResourceId: null,
        queueName: null,
        coreCount: null,
        nodeCount: null,
        wallTimeLimit: null,
        totalPhysicalMemory: null,
        ...this.inputState,
      };
    },
    epolyscatApplicationModuleId() {
      return this.applicationModuleId;
    },
    filteredSaveTargetOptions() {
      const query = this.saveTargetSearch.trim().toLowerCase();
      if (!query) return this.saveTargetOptions;

      return this.saveTargetOptions.filter(option => option.toLowerCase().indexOf(query) !== -1);
    },
    selectedSaveTargetLabel() {
      return this.saveTarget || "Select";
    },
    filteredAllocationOptions() {
      return this.filterResourceOptions(this.allocationOptions, this.allocationSearch);
    },
    filteredComputeResourceOptions() {
      return this.filterResourceOptions(this.computeResourceOptions, this.computeResourceSearch);
    },
    filteredQueueOptions() {
      return this.filterResourceOptions(this.queueOptions, this.queueSearch);
    },
    computeResourceEmptyMessage() {
      if (!this.groupResourceProfileId) return "Select allocation first";
      if (this.resourceReadiness.status === "checking") return "Checking deployments";
      if (this.resourceReadiness.status === "load-error") return "Unable to load deployments";
      if (this.resourceReadiness.status === "no-deployment") return "No deployment available";
      return "Select an available option";
    },
    computeResourceValidationMessage() {
      return this.resourceReadiness.ready
        ? "Compute resource is required"
        : this.resourceReadiness.message;
    },
    queueDisplayName() {
      if (this.queueName) return this.queueName;
      return this.computeResourceId ? "Select queue" : "Select compute resource";
    },
    selectedApplicationDeployment() {
      return this.applicationDeployments.find(
          deployment => deployment.computeHostId === this.computeResourceId
      ) || null;
    },
    selectedComputeResourceLabel() {
      return this.resourceOptionLabel(
          this.computeResourceOptions,
          this.computeResourceId
      );
    },
    deploymentScopeLoaded() {
      return this.loadedApplicationModuleId === this.epolyscatApplicationModuleId
        && this.loadedGroupResourceProfileId === this.groupResourceProfileId
        && this.loadedDeploymentExecutionKind === this.deploymentExecutionKind;
    },
    resourceReadiness() {
      return buildComputeResourceReadiness({
        applicationModuleId: this.epolyscatApplicationModuleId,
        applicationLabel: this.applicationLabel,
        groupResourceProfileId: this.groupResourceProfileId,
        computeResourceId: this.computeResourceId,
        computeResourceLabel: this.selectedComputeResourceLabel,
        eligibleDeployments: this.applicationDeployments,
        scopeLoaded: this.deploymentScopeLoaded,
        loading: this.computeResourcesLoading,
        error: this.computeResourceLoadError,
      });
    },
    resourceReadinessIcon() {
      if (this.resourceReadiness.status === "ready") return "check-circle";
      if (this.resourceReadiness.status === "checking") return "clock";
      if (["selection-required", "allocation-required"].includes(this.resourceReadiness.status)) {
        return "info-circle";
      }
      return "exclamation-circle";
    },
    selectedGroupResourceProfile() {
      return this.groupResourceProfiles.find(
          groupResourceProfile => groupResourceProfile.groupResourceProfileId === this.groupResourceProfileId
      ) || null;
    },
    selectedComputeResourcePolicy() {
      if (!this.selectedGroupResourceProfile || !this.computeResourceId) return null;

      if (this.selectedGroupResourceProfile.getComputeResourcePolicy) {
        return this.selectedGroupResourceProfile.getComputeResourcePolicy(this.computeResourceId);
      }

      return (this.selectedGroupResourceProfile.computeResourcePolicies || [])
          .find(policy => policy.computeResourceId === this.computeResourceId) || null;
    },
    selectedBatchQueueResourcePolicies() {
      if (!this.selectedGroupResourceProfile || !this.computeResourceId) return [];

      if (this.selectedGroupResourceProfile.getBatchQueueResourcePolicies) {
        return this.selectedGroupResourceProfile.getBatchQueueResourcePolicies(this.computeResourceId);
      }

      return (this.selectedGroupResourceProfile.batchQueueResourcePolicies || [])
          .filter(policy => policy.computeResourceId === this.computeResourceId);
    },
    selectedQueueDefault() {
      const option = this.queueOptions.find(queueOption => queueOption.value === this.queueName);
      return option ? option.queueDefault : null;
    },
    defaultQueue() {
      return this.queueOptions.length > 0 ? this.queueOptions[0] : null;
    },
    maxCoreCount() {
      if (!this.selectedQueueDefault) return 0;

      const batchQueueResourcePolicy = this.getBatchQueueResourcePolicy(this.selectedQueueDefault.queueName);
      const queueMax = this.selectedQueueDefault.maxProcessors || Number.MAX_SAFE_INTEGER;
      if (batchQueueResourcePolicy && batchQueueResourcePolicy.maxAllowedCores) {
        return Math.min(batchQueueResourcePolicy.maxAllowedCores, queueMax);
      }
      return queueMax;
    },
    maxNodeCount() {
      if (!this.selectedQueueDefault) return 0;

      const batchQueueResourcePolicy = this.getBatchQueueResourcePolicy(this.selectedQueueDefault.queueName);
      const queueMax = this.selectedQueueDefault.maxNodes || Number.MAX_SAFE_INTEGER;
      if (batchQueueResourcePolicy && batchQueueResourcePolicy.maxAllowedNodes) {
        return Math.min(batchQueueResourcePolicy.maxAllowedNodes, queueMax);
      }
      return queueMax;
    },
    maxWallTimeLimit() {
      if (!this.selectedQueueDefault) return 0;

      const batchQueueResourcePolicy = this.getBatchQueueResourcePolicy(this.selectedQueueDefault.queueName);
      const queueMax = this.selectedQueueDefault.maxRunTime || Number.MAX_SAFE_INTEGER;
      if (batchQueueResourcePolicy && batchQueueResourcePolicy.maxAllowedWalltime) {
        return Math.min(batchQueueResourcePolicy.maxAllowedWalltime, queueMax);
      }
      return queueMax;
    },
  },
  methods: {
    updateRunResources(resources) {
      if (this.disabled) return;
      this.$emit("updateResources", resources);
    },
    selectSaveTarget(option) {
      this.saveTargetSearch = "";
      this.showSaveTargetMenu = false;
      this.$emit("saveTargetSelected", option);
    },
    openResourceMenu(kind) {
      if (this.disabled) return;

      this.loadResourceOptions();
      this.showQueueMenu = false;

      if (kind === "allocation") {
        this.showAllocationMenu = true;
        this.showComputeResourceMenu = false;
      } else if (kind === "compute") {
        this.showComputeResourceMenu = true;
        this.showAllocationMenu = false;
        if (this.computeResourceLoadError && !this.computeResourcesLoading) {
          this.loadComputeResourceOptions();
        }
      }

      this.showSaveTargetMenu = false;
    },
    toggleResourceMenu(kind) {
      if (this.disabled) return;

      if (kind === "allocation") {
        const shouldOpen = !this.showAllocationMenu;
        this.showAllocationMenu = false;
        this.showComputeResourceMenu = false;
        if (shouldOpen) this.openResourceMenu("allocation");
      } else if (kind === "compute") {
        const shouldOpen = !this.showComputeResourceMenu;
        this.showComputeResourceMenu = false;
        this.showAllocationMenu = false;
        if (shouldOpen) this.openResourceMenu("compute");
      }
    },
    onAllocationSearchInput() {
      this.showAllocationMenu = true;
      this.showComputeResourceMenu = false;

      if (this.groupResourceProfileId) {
        this.updateRunResources({
          groupResourceProfileId: null,
          computeResourceId: null,
        });
        this.computeResourceSearch = "";
        this.clearQueueOptions();
        this.loadComputeResourceOptions();
      }
    },
    onComputeResourceSearchInput() {
      this.showComputeResourceMenu = true;
      this.showAllocationMenu = false;

      if (this.computeResourceId) {
        this.computeResourceId = null;
        this.clearQueueOptions();
      }
    },
    async loadResourceOptions() {
      if (this.resourceOptionsLoading || this.resourceOptionsLoaded) return;

      this.resourceOptionsLoading = true;
      this.resourceOptionsError = null;

      try {
        await this.loadAllocationOptions();
        await this.loadComputeResourceOptions();
        this.resourceOptionsLoaded = true;
        this.syncResourceSearches();
      } catch (error) {
        this.resourceOptionsError = error;
      } finally {
        this.resourceOptionsLoading = false;
      }
    },
    async loadAllocationOptions() {
      const groupResourceProfiles = await services.GroupResourceProfileService.list();
      this.groupResourceProfiles = groupResourceProfiles || [];
      this.allocationOptions = this.normalizeGroupResourceProfileOptions(groupResourceProfiles);
      this.syncAllocationSearchToSelection();
    },
    async loadComputeResourceOptions() {
      const requestId = ++this.computeResourceRequestId;
      const applicationModuleId = this.epolyscatApplicationModuleId;
      const groupResourceProfileId = this.groupResourceProfileId;
      const deploymentExecutionKind = this.deploymentExecutionKind;

      this.computeResourceLoadError = null;
      this.loadedApplicationModuleId = null;
      this.loadedGroupResourceProfileId = null;
      this.loadedDeploymentExecutionKind = null;

      if (!groupResourceProfileId) {
        this.computeResourcesLoading = false;
        this.clearComputeResourceOptions();
        return;
      }

      if (!applicationModuleId) {
        this.computeResourcesLoading = false;
        this.clearComputeResourceOptions();
        return;
      }

      this.computeResourcesLoading = true;
      try {
        const [computeResourceNames, applicationDeployments] = await Promise.all([
          services.ComputeResourceService.names(),
          services.ApplicationDeploymentService.list(
              {
                appModuleId: applicationModuleId,
                groupResourceProfileId,
              },
              {ignoreErrors: true}
          ),
        ]);
        if (requestId !== this.computeResourceRequestId) return;

        this.applicationDeployments = filterEligibleApplicationDeployments(
            applicationDeployments,
            deploymentExecutionKind
        );
        const deploymentIds = this.applicationDeployments
            .map(deployment => deployment.computeHostId)
            .filter(Boolean);
        const includedComputeResourceIds = new Set(deploymentIds);

        this.computeResourceOptions = this.normalizeComputeResourceOptions(
            computeResourceNames,
            includedComputeResourceIds
        );
        this.loadedApplicationModuleId = applicationModuleId;
        this.loadedGroupResourceProfileId = groupResourceProfileId;
        this.loadedDeploymentExecutionKind = deploymentExecutionKind;

        if (
            this.computeResourceId
            && !this.computeResourceOptions.find(option => option.value === this.computeResourceId)
        ) {
          this.computeResourceId = null;
          this.computeResourceSearch = "";
        }
        this.syncComputeResourceSearchToSelection();
        if (this.computeResourceId) {
          await this.loadQueueOptionsForSelection({ preserveExisting: true });
        } else {
          this.clearQueueOptions();
        }
      } catch (error) {
        if (requestId !== this.computeResourceRequestId) return;
        this.computeResourceLoadError = error;
        this.clearComputeResourceOptions();
        this.loadedApplicationModuleId = applicationModuleId;
        this.loadedGroupResourceProfileId = groupResourceProfileId;
        this.loadedDeploymentExecutionKind = deploymentExecutionKind;
      } finally {
        if (requestId === this.computeResourceRequestId) {
          this.computeResourcesLoading = false;
        }
      }
    },
    clearComputeResourceOptions() {
      this.applicationDeployments = [];
      this.computeResourceOptions = [];
      if (this.computeResourceId) {
        this.computeResourceId = null;
        this.computeResourceSearch = "";
      }
      this.syncComputeResourceSearchToSelection();
      this.clearQueueOptions();
    },
    async loadQueueOptionsForSelection({ preserveExisting = false } = {}) {
      const applicationDeployment = this.selectedApplicationDeployment;
      if (!applicationDeployment || !applicationDeployment.appDeploymentId) {
        this.clearQueueOptions();
        return;
      }

      this.queueOptionsLoading = true;
      try {
        const queueDefaults = await services.ApplicationDeploymentService.getQueues({
          lookup: applicationDeployment.appDeploymentId,
        }).catch(() => []);
        this.queueOptions = this.normalizeQueueOptions(queueDefaults);

        const existingQueueOption = this.queueOptions.find(option => option.value === this.queueName);
        if (preserveExisting && existingQueueOption) {
          this.applyQueueLimits();
        } else {
          this.setDefaultQueue();
        }
      } finally {
        this.queueOptionsLoading = false;
      }
    },
    clearQueueOptions() {
      this.queueOptions = [];
      this.clearQueueSelection();
    },
    clearQueueSelection() {
      this.queueSearch = "";
      this.showQueueMenu = false;
      this.updateRunResources({
        queueName: null,
        coreCount: null,
        nodeCount: null,
        wallTimeLimit: null,
        totalPhysicalMemory: 0,
      });
    },
    normalizeGroupResourceProfileOptions(groupResourceProfiles) {
      const options = (groupResourceProfiles || []).map(groupResourceProfile => {
        return {
          value: groupResourceProfile.groupResourceProfileId,
          label: groupResourceProfile.groupResourceProfileName || groupResourceProfile.groupResourceProfileId,
        };
      });
      return this.sortResourceOptions(options);
    },
    normalizeComputeResourceOptions(computeResourceNames, includedComputeResourceIds) {
      const names = computeResourceNames || {};
      const options = Object.keys(names)
          .filter(computeResourceId => !includedComputeResourceIds || includedComputeResourceIds.has(computeResourceId))
          .map(computeResourceId => {
            const computeResourceName = names[computeResourceId];
            return {
              value: computeResourceId,
              label: this.getComputeResourceDisplayName(computeResourceName, computeResourceId),
            };
          });
      return this.sortResourceOptions(options);
    },
    normalizeQueueOptions(queueDefaults) {
      return (queueDefaults || [])
          .filter(queueDefault => queueDefault.queueName && this.isQueueAllowedForSelection(queueDefault.queueName))
          .sort((a, b) => {
            if (a.isDefaultQueue) return -1;
            if (b.isDefaultQueue) return 1;
            return a.queueName.localeCompare(b.queueName);
          })
          .map(queueDefault => {
            return {
              value: queueDefault.queueName,
              label: queueDefault.queueName,
              queueDefault,
            };
          });
    },
    isQueueAllowedForSelection(queueName) {
      const computeResourcePolicy = this.selectedComputeResourcePolicy;
      if (
          !computeResourcePolicy
          || !computeResourcePolicy.allowedBatchQueues
          || computeResourcePolicy.allowedBatchQueues.length === 0
      ) {
        return true;
      }
      return computeResourcePolicy.allowedBatchQueues.includes(queueName);
    },
    getBatchQueueResourcePolicy(queueName) {
      return this.selectedBatchQueueResourcePolicies.find(policy => policy.queuename === queueName) || null;
    },
    getComputeResourceDisplayName(computeResourceName, computeResourceId) {
      if (computeResourceName && typeof computeResourceName === "object") {
        return computeResourceName.hostName
            || computeResourceName.name
            || computeResourceName.value
            || computeResourceId;
      }

      return computeResourceName || computeResourceId;
    },
    sortResourceOptions(options) {
      const seen = new Set();
      return options
          .filter(option => {
            if (!option.value || seen.has(option.value)) return false;
            seen.add(option.value);
            return true;
          })
          .sort((a, b) => a.label.localeCompare(b.label));
    },
    filterResourceOptions(options, search) {
      const query = (search || "").trim().toLowerCase();
      if (!query) return options;

      return options.filter(option => {
        return option.label.toLowerCase().indexOf(query) !== -1
            || option.value.toLowerCase().indexOf(query) !== -1;
      });
    },
    resourceOptionLabel(options, value) {
      const option = options.find(candidate => candidate.value === value);
      return option ? option.label : value;
    },
    syncResourceSearches() {
      this.syncAllocationSearchToSelection();
      this.syncComputeResourceSearchToSelection();
    },
    syncAllocationSearchToSelection() {
      if (this.groupResourceProfileId) {
        this.allocationSearch = this.resourceOptionLabel(this.allocationOptions, this.groupResourceProfileId);
      } else if (!this.showAllocationMenu) {
        this.allocationSearch = "";
      }
    },
    syncComputeResourceSearchToSelection() {
      if (this.computeResourceId) {
        this.computeResourceSearch = this.resourceOptionLabel(this.computeResourceOptions, this.computeResourceId);
      } else if (!this.showComputeResourceMenu) {
        this.computeResourceSearch = "";
      }
    },
    selectAllocationOption(option) {
      const previousAllocationId = this.groupResourceProfileId;
      this.groupResourceProfileId = option.value;
      this.allocationSearch = option.label;
      this.showAllocationMenu = false;

      if (previousAllocationId !== option.value) {
        this.computeResourceSearch = "";
        this.updateRunResources({ computeResourceId: null });
        this.clearQueueOptions();
        this.loadComputeResourceOptions();
      }
    },
    selectComputeResourceOption(option) {
      this.computeResourceId = option.value;
      this.computeResourceSearch = option.label;
      this.showComputeResourceMenu = false;
      this.loadQueueOptionsForSelection();
    },
    setDefaultQueue() {
      if (this.defaultQueue) {
        this.selectQueueOption(this.defaultQueue);
      } else {
        this.clearQueueSelection();
      }
    },
    selectQueueOption(option) {
      if (!option) {
        this.clearQueueOptions();
        return;
      }

      this.queueName = option.value;
      this.queueSearch = "";
      this.showQueueMenu = false;
      this.applyQueueDefaults(option.queueDefault);
    },
    toggleQueueMenu() {
      if (this.disabled || this.queueOptions.length === 0) return;

      this.showQueueMenu = !this.showQueueMenu;
      this.showAllocationMenu = false;
      this.showComputeResourceMenu = false;
      this.showSaveTargetMenu = false;
    },
    applyQueueDefaults(queueDefault) {
      if (!queueDefault) return;

      this.updateRunResources({
        coreCount: this.getDefaultCoreCount(queueDefault),
        nodeCount: this.getDefaultNodeCount(queueDefault),
        wallTimeLimit: this.getDefaultWallTimeLimit(queueDefault),
        totalPhysicalMemory: 0,
      });
      this.applyQueueLimits();
    },
    getDefaultCoreCount(queueDefault) {
      return this.applyQueuePolicyLimit(
          queueDefault.defaultCPUCount || 1,
          queueDefault,
          "maxAllowedCores",
          "maxProcessors"
      );
    },
    getDefaultNodeCount(queueDefault) {
      return this.applyQueuePolicyLimit(
          queueDefault.defaultNodeCount || 1,
          queueDefault,
          "maxAllowedNodes",
          "maxNodes"
      );
    },
    getDefaultWallTimeLimit(queueDefault) {
      return this.applyQueuePolicyLimit(
          queueDefault.defaultWalltime || 1,
          queueDefault,
          "maxAllowedWalltime",
          "maxRunTime"
      );
    },
    applyQueuePolicyLimit(value, queueDefault, policyField, queueField) {
      const batchQueueResourcePolicy = this.getBatchQueueResourcePolicy(queueDefault.queueName);
      const queueLimit = queueDefault[queueField] || Number.MAX_SAFE_INTEGER;
      const policyLimit = batchQueueResourcePolicy && batchQueueResourcePolicy[policyField]
          ? batchQueueResourcePolicy[policyField]
          : Number.MAX_SAFE_INTEGER;
      return Math.min(value, queueLimit, policyLimit);
    },
    applyQueueLimits() {
      if (!this.selectedQueueDefault) return;

      const resources = {};
      if (this.maxCoreCount > 0 && this.coreCount > this.maxCoreCount) {
        resources.coreCount = this.maxCoreCount;
      }
      if (this.maxNodeCount > 0 && this.nodeCount > this.maxNodeCount) {
        resources.nodeCount = this.maxNodeCount;
      }
      if (this.maxWallTimeLimit > 0 && this.wallTimeLimit > this.maxWallTimeLimit) {
        resources.wallTimeLimit = this.maxWallTimeLimit;
      }
      if (Object.keys(resources).length > 0) {
        this.updateRunResources(resources);
      }
    },
    toggleQueueResourceLock() {
      if (this.disabled) return;

      this.queueResourceLockEnabled = !this.queueResourceLockEnabled;
      if (this.queueResourceLockEnabled) {
        this.onNodeCountInput();
      }
    },
    onNodeCountInput() {
      if (
          this.queueResourceLockEnabled
          && this.selectedQueueDefault
          && this.selectedQueueDefault.cpuPerNode > 0
      ) {
        const nodeCount = parseInt(this.nodeCount);
        if (nodeCount > 0) {
          this.coreCount = Math.min(
              nodeCount * this.selectedQueueDefault.cpuPerNode,
              this.maxCoreCount || Number.MAX_SAFE_INTEGER
          );
        }
      }
    },
    onCoreCountInput() {
      if (
          this.queueResourceLockEnabled
          && this.selectedQueueDefault
          && this.selectedQueueDefault.cpuPerNode > 0
      ) {
        const coreCount = parseInt(this.coreCount);
        if (coreCount > 0) {
          this.nodeCount = Math.min(
              Math.ceil(coreCount / this.selectedQueueDefault.cpuPerNode),
              this.maxNodeCount || Number.MAX_SAFE_INTEGER
          );
        }
      }
    },
    async refreshResourceSelection() {
      if (!this.resourceOptionsLoaded) {
        await this.loadResourceOptions();
      } else {
        this.syncResourceSearches();
        await this.loadComputeResourceOptions();
      }
    },
  },
  watch: {
    run(newRun, oldRun) {
      if (newRun !== oldRun) {
        this.refreshResourceSelection();
      }
    },
    epolyscatApplicationModuleId() {
      this.loadComputeResourceOptions();
    },
    deploymentExecutionKind() {
      this.loadComputeResourceOptions();
    },
    resourceReadiness: {
      handler(readiness) {
        this.$emit("readinessChanged", { ...readiness });
      },
      immediate: true,
      deep: true,
    },
  },
  async mounted() {
    await this.$store.dispatch("settings/fetchSettings");
    await this.loadResourceOptions();
  },
};
</script>

<style scoped>
.resource-settings-grid {
  column-gap: 141px;
  display: grid;
  grid-template-columns: minmax(0, 438px) minmax(0, 438px);
  margin-top: 64px;
  max-width: 1017px;
  row-gap: 14px;
}

.resource-fields {
  display: grid;
  gap: 14px;
}

.resource-row {
  align-items: start;
  display: grid;
  gap: 10px;
  grid-template-columns: 86px 260px;
}

.resource-row label,
.settings-label {
  color: #000000;
  font-size: 16px;
  font-weight: 400;
  line-height: 38px;
}

.save-to-dropdown {
  position: relative;
  width: 208px;
}

.resource-search-dropdown {
  position: relative;
  width: 260px;
}

.resource-search-input {
  border-color: #d6d6d6;
  border-radius: 6px;
  color: #000000;
  font-size: 15px;
  font-weight: 600;
  height: 38px;
  padding-right: 36px;
}

.resource-search-input.is-valid {
  background-position: right 38px center;
  padding-right: 62px;
}

.resource-readiness {
  align-items: flex-start;
  color: #455967;
  display: flex;
  font-size: 12px;
  font-weight: 600;
  gap: 6px;
  line-height: 1.35;
  margin-top: 7px;
}

.resource-readiness .b-icon {
  flex: 0 0 auto;
  margin-top: 2px;
}

.resource-readiness.status-ready {
  color: #21643f;
}

.resource-readiness.status-load-error,
.resource-readiness.status-no-deployment,
.resource-readiness.status-configuration-missing,
.resource-readiness.status-invalid-selection {
  color: #9b2f2f;
}

.resource-search-toggle {
  align-items: center;
  background: transparent;
  border: 0;
  color: #555555;
  display: inline-flex;
  height: 38px;
  justify-content: center;
  position: absolute;
  right: 0;
  top: 0;
  width: 34px;
}

.resource-search-menu {
  background: #ffffff;
  border: 1px solid #d6d6d6;
  box-shadow: 0 4px 10px rgba(0, 0, 0, 0.12);
  left: 0;
  max-height: 260px;
  overflow-y: auto;
  position: absolute;
  top: 38px;
  width: 260px;
  z-index: 14;
}

.resource-menu-search {
  border-radius: 6px;
  margin: 8px;
  width: calc(100% - 16px);
}

.resource-option {
  background: #ffffff;
  border: 0;
  border-top: 1px solid #eeeeee;
  color: #000000;
  display: grid;
  gap: 2px;
  min-height: 40px;
  padding: 8px 16px;
  text-align: left;
  width: 100%;
}

.resource-option:hover,
.resource-option.active {
  background: #f2f6f8;
}

.resource-option-label {
  font-size: 14px;
  font-weight: 700;
  line-height: 1.2;
}

.resource-empty-state {
  color: #666666;
  font-size: 13px;
  padding: 12px 16px;
}

.save-to-trigger {
  align-items: center;
  background: #ffffff;
  border: 1px solid #d6d6d6;
  border-radius: 6px;
  color: #999999;
  display: flex;
  font-size: 16px;
  font-weight: 700;
  height: 38px;
  justify-content: space-between;
  padding: 0 12px 0 78px;
  width: 100%;
}

.save-to-menu {
  background: #ffffff;
  border: 1px solid #d6d6d6;
  box-shadow: 0 4px 10px rgba(0, 0, 0, 0.12);
  left: 0;
  position: absolute;
  top: 38px;
  width: 208px;
  z-index: 8;
}

.save-to-search {
  border-radius: 6px;
  margin: 8px;
  width: calc(100% - 16px);
}

.save-to-option {
  background: #ffffff;
  border: 0;
  border-top: 1px solid #eeeeee;
  color: #000000;
  display: block;
  font-size: 16px;
  min-height: 40px;
  padding: 8px 16px;
  text-align: left;
  width: 100%;
}

.save-to-option:hover {
  background: #f2f6f8;
}

.settings-label {
  grid-column: 1;
  margin-top: 12px;
}

.queue-settings-card {
  background: #f4f4f4;
  border-radius: 6px;
  display: grid;
  gap: 18px;
  grid-column: 1 / span 2;
  margin-left: 86px;
  min-height: 120px;
  padding: 16px 28px 18px;
  position: relative;
  width: calc(100% - 86px);
}

.queue-card-header {
  align-items: start;
  column-gap: 16px;
  display: grid;
  grid-template-columns: minmax(0, 1fr) auto;
}

.queue-selector {
  justify-self: end;
  position: relative;
}

.queue-check {
  align-items: center;
  background: #226597;
  border: 0;
  border-radius: 7px;
  color: #ffffff;
  display: inline-flex;
  font-size: 24px;
  height: 32px;
  justify-content: center;
  width: 32px;
}

.queue-menu {
  background: #ffffff;
  border: 1px solid #d6d6d6;
  box-shadow: 0 4px 10px rgba(0, 0, 0, 0.12);
  position: absolute;
  right: 0;
  top: 38px;
  width: 220px;
  z-index: 14;
}

.queue-search {
  border-radius: 6px;
  margin: 8px;
  width: calc(100% - 16px);
}

.queue-option {
  background: #ffffff;
  border: 0;
  border-top: 1px solid #eeeeee;
  color: #000000;
  display: block;
  font-size: 14px;
  font-weight: 700;
  min-height: 40px;
  padding: 8px 16px;
  text-align: left;
  width: 100%;
}

.queue-option:hover,
.queue-option.active {
  background: #f2f6f8;
}

.queue-empty-state {
  color: #666666;
  font-size: 13px;
  padding: 12px 16px;
}

.queue-check:disabled {
  background: #d0d0d0;
  color: #777777;
}

.queue-name {
  color: #000000;
  font-size: 16px;
  line-height: 32px;
  min-height: 32px;
  overflow-wrap: anywhere;
}

.queue-fields {
  align-items: start;
  column-gap: 28px;
  display: grid;
  grid-template-columns: minmax(92px, 1fr) 32px minmax(92px, 1fr) minmax(118px, 1fr) minmax(176px, 1.35fr);
  row-gap: 14px;
}

.queue-fields label {
  display: grid;
  gap: 7px;
  min-width: 0;
}

.queue-lock-button {
  align-items: center;
  background: #ffffff;
  border: 1px solid #d6d6d6;
  border-radius: 6px;
  color: #555555;
  display: inline-flex;
  font-size: 15px;
  height: 32px;
  justify-content: center;
  margin-top: 0;
  min-width: 32px;
  padding: 0;
  width: 32px;
}

.queue-lock-button .b-icon {
  flex: 0 0 16px;
  height: 16px;
  width: 16px;
}

.queue-lock-button.active {
  background: #226597;
  border-color: #226597;
  color: #ffffff;
}

.queue-lock-button:disabled {
  background: #eeeeee;
  border-color: #dddddd;
  color: #999999;
}

.queue-input {
  background: #ffffff;
  border-color: #d6d6d6;
  border-radius: 6px;
  color: #000000;
  font-size: 14px;
  height: 32px;
  min-width: 0;
  text-align: center;
  width: 100%;
}

.memory-input {
  min-width: 0;
}

.queue-fields span {
  color: #000000;
  font-size: 14px;
  line-height: 1.2;
  min-height: 18px;
  overflow-wrap: break-word;
  white-space: normal;
}

.resource-validation-copy,
.native-resource-selectors {
  height: 1px;
  left: -10000px;
  opacity: 0;
  overflow: hidden;
  pointer-events: none;
  position: absolute;
  top: auto;
  width: 1px;
}

@media (max-width: 920px) {
  .resource-settings-grid {
    grid-template-columns: 1fr;
    row-gap: 20px;
  }

  .queue-settings-card {
    grid-column: 1;
    margin-left: 0;
    width: 100%;
  }

  .queue-fields {
    grid-template-columns: repeat(2, minmax(90px, 1fr));
  }
}
</style>
