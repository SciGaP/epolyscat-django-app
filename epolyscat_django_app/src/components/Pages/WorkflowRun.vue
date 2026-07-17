<template>
  <div class="workflow-run-page">
    <main class="workflow-run-content">
      <b-breadcrumb class="workflow-breadcrumb">
        <b-breadcrumb-item to="/runs">Experiments</b-breadcrumb-item>
        <b-breadcrumb-item active>{{ run.name }}</b-breadcrumb-item>
      </b-breadcrumb>

      <header class="workflow-heading">
        <div>
          <h1>{{ run.name }}</h1>
          <div class="workflow-subtitle">Workflow/{{ activeStageLabel }}</div>
        </div>
        <b-button
            variant="outline-primary"
            class="run-type-switch-button"
            v-on:click="switchToModuleRun"
        >
          Switch to Module Run
        </b-button>
      </header>

      <nav class="workflow-stepper" aria-label="Run workflow">
        <button
            v-for="(stage, index) in workflowStages"
            :key="stage.id"
            type="button"
            class="workflow-step"
            :class="{ active: stage.id === activeStageId, complete: stage.state === 'complete', disabled: stage.disabled }"
            :disabled="stage.disabled"
            :aria-current="stage.id === activeStageId ? 'step' : null"
            v-on:click="selectStage(stage)"
        >
          <span class="workflow-step-index">{{ index + 1 }}</span>
          <span class="workflow-step-label">{{ stage.label }}</span>
        </button>
      </nav>

      <section class="workflow-section">
        <h2>Applications</h2>
        <div class="application-choice-row">
          <button
              v-for="application in applications"
              :key="application.id"
              type="button"
              class="application-choice"
              :class="{ active: workflowApplication === application.id }"
              v-on:click="workflowApplication = application.id"
          >
            {{ application.label }}
          </button>
        </div>
      </section>

      <section class="workflow-section">
        <h2>{{ workflowApplication }} Parameters</h2>
        <div class="workflow-parameter-card">
          <div
              v-for="parameter in activeApplication.parameters"
              :key="parameter.name"
              class="workflow-parameter-row"
          >
            <label>{{ parameter.label }}</label>
            <template v-if="parameter.kind === 'file'">
              <div class="workflow-file-actions">
                <b-button variant="dark" class="workflow-storage-button" v-b-modal.workflow-input-storage>
                  Select file from storage
                </b-button>
                <span class="workflow-or">OR</span>
                <label class="workflow-computer-button">
                  Select file from computer
                  <input type="file" v-on:change="onComputerFileSelected">
                </label>
              </div>
              <div v-if="selectedWorkflowFileName" class="workflow-selected-file">
                {{ selectedWorkflowFileName }}
              </div>
            </template>
            <b-form-select
                v-else
                v-model="workflowMetadata[parameter.name]"
                :options="['Yes', 'No']"
            />
          </div>
        </div>
        <UserStorage
            id="workflow-input-storage"
            selectionType="files"
            :canSelectMultiple="false"
            modalTitle="Select workflow input file"
            :onlyShowEditable="false"
            @filesSelected="onStorageFilesSelected"
        />
      </section>

      <section class="workflow-section">
        <div class="workflow-section-title-row">
          <h2>Input Files</h2>
          <div class="workflow-view-toggle">
            <b-icon icon="list-ul" />
            <b-icon icon="card-list" />
          </div>
        </div>
        <div class="workflow-input-tabs" role="tablist" aria-label="Workflow input files">
          <button type="button" class="workflow-input-tab active">input_data</button>
          <button type="button" class="workflow-input-tab">input_file (.inp)</button>
        </div>
        <pre class="workflow-code-preview">{{ workflowInputPreview }}</pre>
      </section>

      <section class="workflow-section">
        <h2>Resource</h2>
        <RunResourceSettings
            :run="run"
            :input-state="inputState"
            v-on:updateResources="updateResources"
        />
      </section>

      <div class="workflow-actions">
        <button-overlay :show="processing">
          <b-button variant="outline-primary" :disabled="processing" v-on:click="onSave(false)">
            Save
          </b-button>
        </button-overlay>
        <button-overlay :show="processing">
          <b-button variant="primary" :disabled="processing" v-on:click="onSave(true)">
            Submit
          </b-button>
        </button-overlay>
      </div>
    </main>
  </div>
</template>

<script>
import router from "@/router";
import { eventBus } from "@/event-bus";
import ButtonOverlay from "@/components/overlay/button-overlay";
import RunResourceSettings from "@/components/blocks/RunResourceSettings";
import UserStorage from "@/components/overlay/UserStorage";

export default {
  name: "WorkflowRun",
  components: { ButtonOverlay, RunResourceSettings, UserStorage },
  data() {
    return {
      runMode: "workflow",
      activeStageId: "Data_Gen",
      workflowApplication: "OpenMolcas",
      utilityApplication: "CnvMath",
      selectedWorkflowFile: null,
      selectedWorkflowFilePayload: null,
      workflowMetadata: {
        gpu_version: "Yes",
        printing_orbitals: "Yes",
      },
      processing: false,
      run: {
        name: `Workflow Run on ${(new Date(Date.now())).toLocaleString()}`,
        description: "",
        status: "UNSUBMITTED",
        groupResourceProfileId: null,
        computeResourceId: null,
        queueName: null,
        coreCount: null,
        nodeCount: null,
        wallTimeLimit: null,
        totalPhysicalMemory: null,
      },
      workflowStages: [
        { id: "Data_Gen", label: "Data Generation", state: "active" },
        { id: "ePolyScat_Run", label: "ePolyScat Run", state: "pending" },
        { id: "Analysis", label: "Post-processing", state: "pending", disabled: true },
        { id: "Visualization", label: "Visualization", state: "pending", disabled: true, localOnly: true },
      ],
      applications: [
        {
          id: "Gaussian16",
          label: "Gaussian16",
          parameters: [
            { name: "gaussian_input", label: "Input file", kind: "file" },
            { name: "gpu_version", label: "GPU_Version?", kind: "choice" },
          ],
        },
        {
          id: "OpenMolcas",
          label: "OpenMolcas",
          parameters: [
            { name: "openmolcas_input", label: "Input file", kind: "file" },
            { name: "openmolcas_optional", label: "Optional files", kind: "file" },
            { name: "printing_orbitals", label: "Printing of orbitals", kind: "choice" },
          ],
        },
      ],
    };
  },
  computed: {
    activeStageLabel() {
      const stage = this.workflowStages.find(stage => stage.id === this.activeStageId);
      return stage ? stage.label : this.activeStageId;
    },
    activeApplication() {
      return this.applications.find(application => application.id === this.workflowApplication) || this.applications[0];
    },
    activeWorkflowInputName() {
      if (this.activeStageId === "Data_Gen") {
        return this.workflowApplication === "Gaussian16" ? "Gaussian_Inputs" : "Molcas_Inputs";
      }
      if (this.activeStageId === "Analysis" && this.utilityApplication === "MoldenMerge") {
        return "molden.dat";
      }
      return "ePolyScat_Input_Data";
    },
    selectedWorkflowFileName() {
      return this.selectedWorkflowFile ? this.selectedWorkflowFile.name : "";
    },
    workflowInputPreview() {
      return [
        `# ${this.run.name}`,
        "Calculation_Type = WORKFLOW",
        `Application_Workflow = ${this.activeStageId}`,
        this.activeStageId === "Data_Gen" ? `Data_Gen = ${this.workflowApplication}` : null,
        this.activeStageId === "Analysis" ? `Application_Utility = ${this.utilityApplication}` : null,
        `Input_File = ${this.selectedWorkflowFileName || "target"}`,
      ].filter(Boolean).join("\n");
    },
    inputState() {
      return {
        groupResourceProfileId: this.run.groupResourceProfileId === null ? null : !!this.run.groupResourceProfileId,
        computeResourceId: this.run.computeResourceId === null ? null : !!this.run.computeResourceId,
        queueName: this.run.queueName === null ? null : !!this.run.queueName,
        coreCount: this.run.coreCount === null ? null : this.run.coreCount > 0,
        nodeCount: this.run.nodeCount === null ? null : this.run.nodeCount > 0,
        wallTimeLimit: this.run.wallTimeLimit === null ? null : this.run.wallTimeLimit > 0,
        totalPhysicalMemory: this.run.totalPhysicalMemory === null ? null : `${this.run.totalPhysicalMemory}`.length >= 1,
      };
    },
  },
  methods: {
    switchToModuleRun() {
      router.push({
        path: "/runs/new",
        query: this.$route.query
      });
    },
    selectStage(stage) {
      if (!stage.disabled) this.activeStageId = stage.id;
    },
    updateResources(resources) {
      Object.assign(this.run, resources);
    },
    onComputerFileSelected(event) {
      const file = (event.target.files || [])[0];
      if (file) {
        file.isFromComputer = true;
        this.selectedWorkflowFile = file;
        const reader = new FileReader();
        reader.onload = () => {
          this.selectedWorkflowFilePayload = {
            name: file.name,
            contents: reader.result,
            isPlaintext: true,
            deleted: false,
          };
        };
        reader.readAsText(file);
      }
      event.target.value = "";
    },
    onStorageFilesSelected(files) {
      const file = files && files[0];
      if (file) {
        this.selectedWorkflowFile = {
          ...file,
          isFromComputer: false,
        };
        this.selectedWorkflowFilePayload = null;
      }
    },
    buildInputs() {
      const generatedFile = {
        name: "workflow_input.inp",
        contents: this.workflowInputPreview,
        isPlaintext: true,
        deleted: false,
      };
      const files = [generatedFile];
      if (this.selectedWorkflowFilePayload) {
        files.push(this.selectedWorkflowFilePayload);
      }
      if (this.selectedWorkflowFile && this.selectedWorkflowFile.dataProductURI) {
        files.push({
          name: this.selectedWorkflowFile.name,
          dataProductURI: this.selectedWorkflowFile.dataProductURI,
          deleted: false,
        });
      }

      const inputs = [
        { type: "radio input", name: "Calculation_Type", value: "WORKFLOW" },
        { type: "radio input", name: "Application_Workflow", value: this.activeStageId },
      ];
      if (this.activeStageId === "Data_Gen") {
        inputs.push({ type: "radio input", name: "Data_Gen", value: this.workflowApplication });
      }
      if (this.activeStageId === "Analysis") {
        inputs.push({ type: "radio input", name: "Application_Utility", value: this.utilityApplication });
      }
      inputs.push({ type: "files", name: this.activeWorkflowInputName, files });
      return inputs;
    },
    async onSave(submit = false) {
      this.processing = true;
      try {
        let run = await this.$store.dispatch("run/createRun", {
          ...this.run,
          runMode: "workflow",
          workflowStage: this.activeStageId,
          workflowApplication: this.activeStageId === "Data_Gen" ? this.workflowApplication : "",
          utilityApplication: this.activeStageId === "Analysis" ? this.utilityApplication : "",
          workflowMetadata: {
            ...this.workflowMetadata,
            selected_file: this.selectedWorkflowFileName || "target",
            code: this.workflowInputPreview,
          },
          inputs: this.buildInputs(),
          viewIds: [],
        });

        if (submit) {
          run = await this.$store.dispatch("run/submitRun", { runId: run.id });
        }

        router.push(`/runs/${run.id}`);
      } catch (error) {
        eventBus.$emit("error", { name: `Error while trying to ${submit ? "submit" : "save"} workflow run`, error });
      } finally {
        this.processing = false;
      }
    },
  },
};
</script>

<style scoped>
.workflow-run-page {
  background: #ffffff;
  color: #000000;
  min-height: 100%;
}

.workflow-run-content {
  max-width: 1120px;
  padding: 18px 32px 90px 44px;
}

.workflow-breadcrumb {
  font-size: 16px;
  margin-bottom: 18px;
}

.workflow-heading {
  align-items: flex-start;
  display: flex;
  justify-content: space-between;
  margin-bottom: 26px;
}

.workflow-heading h1 {
  font-size: 28px;
  font-weight: 700;
  line-height: 1.2;
  margin: 0 0 6px;
}

.workflow-subtitle {
  color: #555555;
  font-size: 17px;
  font-weight: 700;
}

.run-type-switch-button {
  border-radius: 6px;
  font-weight: 700;
}

.workflow-stepper {
  align-items: center;
  display: flex;
  justify-content: center;
  margin-bottom: 34px;
}

.workflow-step {
  align-items: center;
  background: transparent;
  border: 0;
  color: #555555;
  display: flex;
  font-size: 14px;
  font-weight: 700;
  gap: 8px;
  padding: 0 18px;
  position: relative;
}

.workflow-step + .workflow-step::before {
  border-top: 2px dotted #b8b8b8;
  content: "";
  left: -18px;
  position: absolute;
  top: 16px;
  width: 36px;
}

.workflow-step-index {
  align-items: center;
  background: #ffffff;
  border: 2px solid #b8b8b8;
  border-radius: 50%;
  display: inline-flex;
  height: 32px;
  justify-content: center;
  width: 32px;
}

.workflow-step.complete .workflow-step-index {
  background: #2f9e44;
  border-color: #2f9e44;
  color: #ffffff;
}

.workflow-step.active .workflow-step-index {
  background: #226597;
  border-color: #226597;
  color: #ffffff;
}

.workflow-step.disabled {
  cursor: not-allowed;
  opacity: 0.55;
}

.workflow-section {
  margin-bottom: 34px;
}

.workflow-section h2,
.workflow-section-title-row h2 {
  font-size: 17px;
  font-weight: 700;
  line-height: 1.3;
  margin: 0 0 12px;
}

.workflow-section-title-row {
  align-items: center;
  display: flex;
  justify-content: space-between;
}

.workflow-view-toggle {
  color: #555555;
  display: flex;
  gap: 12px;
}

.application-choice-row {
  display: grid;
  gap: 14px;
  max-width: 720px;
}

.application-choice {
  background: #eeeeee;
  border: 0;
  border-radius: 6px;
  color: #000000;
  font-size: 17px;
  font-weight: 700;
  min-height: 48px;
  text-align: left;
  padding: 0 22px;
}

.application-choice.active {
  background: #226597;
  color: #ffffff;
}

.workflow-parameter-card {
  border: 1px solid #dddddd;
  border-radius: 6px;
  max-width: 760px;
  padding: 18px;
}

.workflow-parameter-row {
  display: grid;
  gap: 12px;
  grid-template-columns: 170px minmax(0, 1fr);
  margin-bottom: 16px;
}

.workflow-parameter-row:last-child {
  margin-bottom: 0;
}

.workflow-parameter-row label {
  font-weight: 700;
  margin: 0;
}

.workflow-file-actions {
  align-items: center;
  display: grid;
  gap: 12px;
  grid-template-columns: minmax(180px, max-content) 36px minmax(220px, 1fr);
}

.workflow-storage-button {
  border-radius: 6px;
}

.workflow-or {
  color: #555555;
  font-size: 12px;
  font-weight: 700;
  text-align: center;
}

.workflow-computer-button {
  align-items: center;
  border: 1px dashed #aaaaaa;
  border-radius: 6px;
  color: #000000;
  cursor: pointer;
  display: flex;
  font-weight: 700;
  justify-content: center;
  margin: 0;
  min-height: 38px;
}

.workflow-computer-button input {
  height: 1px;
  opacity: 0;
  overflow: hidden;
  position: absolute;
  width: 1px;
}

.workflow-selected-file {
  color: #555555;
  font-size: 13px;
  grid-column: 2;
}

.workflow-input-tabs {
  border: 1px solid #dddddd;
  border-radius: 6px;
  display: inline-flex;
  gap: 8px;
  padding: 16px;
}

.workflow-input-tab {
  background: #eeeeee;
  border: 0;
  border-radius: 6px;
  color: #000000;
  font-size: 16px;
  font-weight: 700;
  min-height: 36px;
  padding: 8px 18px;
}

.workflow-input-tab.active {
  background: #226597;
  color: #ffffff;
}

.workflow-code-preview {
  background: #f7f7f7;
  border: 1px solid #dddddd;
  border-radius: 6px;
  color: #1d1d1d;
  font-size: 13px;
  line-height: 1.55;
  margin: 18px 0 0;
  max-width: 760px;
  min-height: 150px;
  padding: 16px;
  white-space: pre-wrap;
}

.workflow-actions {
  display: flex;
  gap: 12px;
  justify-content: flex-end;
  margin-top: 46px;
  max-width: 760px;
}

@media (max-width: 840px) {
  .workflow-run-content {
    padding: 16px 18px 70px;
  }

  .workflow-heading,
  .workflow-section-title-row {
    display: block;
  }

  .workflow-stepper,
  .workflow-step {
    align-items: flex-start;
    flex-direction: column;
  }

  .workflow-step + .workflow-step::before {
    display: none;
  }

  .workflow-parameter-row,
  .workflow-file-actions {
    grid-template-columns: 1fr;
  }

  .workflow-selected-file {
    grid-column: 1;
  }
}
</style>
