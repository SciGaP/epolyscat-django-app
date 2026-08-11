<template>
  <div class="view-run-page">
    <main class="view-run-content">
      <b-breadcrumb class="view-run-breadcrumb">
        <b-breadcrumb-item to="/runs">Experiments</b-breadcrumb-item>
        <b-breadcrumb-item active>{{ run ? run.name : "Run" }}</b-breadcrumb-item>
      </b-breadcrumb>

      <header class="view-run-heading">
        <div>
          <h1>{{ run ? run.name : "Run" }}</h1>
          <div class="view-run-subtitle">
            <span v-if="presentation.mode === 'workflow'">{{ workflowSubtitle }}</span>
            <span v-else>Modules/EPOLYSCAT_DMAT</span>
          </div>
        </div>
        <div class="view-run-heading-actions">
          <div class="view-run-status" v-if="statusBadges.length">
            <span
                v-for="badge in statusBadges"
                :key="badge.label"
                :title="badge.title || null"
            >
              <strong>{{ badge.label }}</strong>
              {{ badge.value }}
            </span>
          </div>
          <b-button
              v-if="workflowContinuation.eligible"
              variant="primary"
              :disabled="workflowContinuationLoading"
              v-on:click="continueInWorkflow"
          >
            {{ workflowContinuationLoading ? "Creating Workflow..." : "Continue in Workflow" }}
          </b-button>
        </div>
      </header>

      <section
          v-if="workflowContinuation.eligible && workflowContinuation.nextStagePreview"
          class="workflow-continuation-preview"
          aria-label="Workflow continuation preview"
      >
        <div v-for="item in continuationPreviewItems" :key="item.label">
          <span class="workflow-status-label">{{ item.label }}</span>
          <strong>{{ item.value }}</strong>
          <span v-if="item.detail" class="workflow-status-application">
            {{ item.detail }}
          </span>
        </div>
      </section>

      <nav
          v-if="presentation.mode === 'workflow'"
          class="workflow-stepper"
          aria-label="Run workflow"
      >
        <span
            v-for="(stage, index) in presentation.stages"
            :key="stage.id"
            class="workflow-step"
            :class="stage.state"
        >
          <span class="workflow-step-index">{{ index + 1 }}</span>
          <span class="workflow-step-copy">
            <router-link
                v-if="workflowStageRoute(stage, index)"
                class="workflow-step-label workflow-step-link"
                :class="{ 'workflow-step-configure-link': workflowStageIsConfigurable(stage, index) }"
                :to="workflowStageRoute(stage, index)"
            >
              {{ stage.label }}
            </router-link>
            <span v-else class="workflow-step-label">{{ stage.label }}</span>
            <span v-if="workflowStageStatusText(stage)" class="workflow-step-state">
              {{ workflowStageStatusText(stage) }}
            </span>
            <span v-if="stage.application" class="workflow-step-application">
              {{ stage.application }}
            </span>
          </span>
        </span>
      </nav>

      <section v-if="presentation.mode === 'workflow'" class="workflow-status-panel">
        <div v-for="item in workflowPanelItems" :key="item.label">
          <span class="workflow-status-label">{{ item.label }}</span>
          <strong>{{ item.value }}</strong>
          <span v-if="item.detail" class="workflow-status-application">
            {{ item.detail }}
          </span>
        </div>
      </section>

      <div class="view-run-layout" :class="{ 'with-plot-panel': showPlotPanel }">
        <section class="view-run-main">
          <section class="target-states-matrix" v-if="showTargetStates">
            <div class="target-column">
              <h2>Input</h2>
              <div v-if="targetInputColumns.length" class="target-input-grid">
                <div
                    v-for="(column, columnIndex) in targetInputColumns"
                    :key="`input-${columnIndex}`"
                    class="target-input-column"
                >
                  <button
                      v-for="file in column"
                      :key="file"
                      type="button"
                      class="file-chip"
                      :class="{ active: selectedFile === file }"
                      v-on:click="selectedFile = file"
                  >
                    {{ file }}
                  </button>
                </div>
              </div>
              <div v-else class="target-empty-state">No input files required</div>
            </div>
            <div class="target-column">
              <h2>Output</h2>
              <div v-if="targetOutputFiles.length" class="target-output-list">
                <button
                    v-for="file in targetOutputFiles"
                    :key="file"
                    type="button"
                    class="file-chip"
                    v-on:click="selectedFile = file"
                >
                  {{ file }}
                </button>
              </div>
              <div v-else class="target-empty-state">No output files available yet</div>
            </div>
          </section>

          <section class="module-parameters" v-if="presentation.mode !== 'workflow'">
            <div class="module-tabs">
              <button type="button" class="module-tab active">Title</button>
              <button type="button" class="module-tab">Settings</button>
              <button type="button" class="module-tab">Target states</button>
              <button type="button" class="module-tab">Partial Waves</button>
            </div>
            <div class="module-parameter-grid">
              <div
                  v-for="parameter in presentation.parameters"
                  :key="parameter.label"
                  class="module-parameter-row"
              >
                <span>{{ parameter.label }}</span>
                <strong>{{ parameter.value }}</strong>
              </div>
            </div>
          </section>

          <section class="run-code-viewer">
            <div class="code-viewer-toolbar">
              <span>{{ selectedFile }}</span>
              <div class="file-control-row">
                <label for="run-file-selector">File</label>
                <select
                    id="run-file-selector"
                    v-model="selectedFile"
                    class="file-selector"
                >
                  <optgroup
                      v-for="group in resolvedFileGroups"
                      :key="group.label"
                      :label="group.label"
                  >
                    <option v-for="file in group.files" :key="`${group.label}-${file}`" :value="file">
                      {{ file }}
                    </option>
                  </optgroup>
                </select>
              </div>
            </div>
            <pre>{{ codePreview }}</pre>
          </section>

          <section class="view-run-resource">
            <h2>Resource</h2>
            <RunResource v-if="run" :run="run" :viewing="true" />
          </section>
        </section>

        <aside class="plot-panel" v-if="showPlotPanel">
          <div class="plot-canvas">
            <img v-if="plotImageUrl" :src="plotImageUrl" alt="Plot result">
            <div v-else class="plot-placeholder">
              <span v-if="plotLoading">Creating plot...</span>
              <span v-else-if="plotError">{{ plotError }}</span>
              <span v-else-if="!hasPlottableOutputFiles">No plottable output file is available yet.</span>
              <span v-else>Select an output file and create a plot.</span>
            </div>
          </div>
          <label>
            File
            <b-form-select v-model="plotForm.file" :options="plotFileOptions" />
          </label>
          <label>
            x-axis
            <b-form-input v-model="plotForm.x_axis" />
          </label>
          <label>
            y-axis
            <b-form-input v-model="plotForm.y_axis" />
          </label>
          <label>
            flags
            <b-form-input v-model="plotForm.flags" />
          </label>
          <b-button
              variant="primary"
              :disabled="!canCreatePlot || plotLoading"
              v-on:click="createPlot"
          >
            {{ plotLoading ? "Plotting" : "Plot" }}
          </b-button>
          <pre v-if="plotOutput" class="plot-output">{{ plotOutput }}</pre>
          <div v-if="plotUserGuidance" class="plot-guidance">
            {{ plotUserGuidance }}
          </div>
        </aside>
      </div>
    </main>
  </div>
</template>

<script>
import { eventBus } from "@/event-bus";
import { InputService, PlotService, RunService } from "@/service/epolyscat-service";
import RunResource from "@/components/blocks/RunResource";

export default {
  name: "ViewRun",
  components: { RunResource },
  data() {
    return {
      run: null,
      presentation: {
        mode: "module",
        subtitle: "Modules/EPOLYSCAT_DMAT",
        file_groups: [],
        target_states: { input_columns: [], output_files: [] },
        parameters: [],
        plot: { file: "", x_axis: "0", y_axis: "1", flags: "-linY" },
        plottable_file_names: [],
        stages: [],
        code: "",
      },
      selectedFile: "target",
      viewables: [],
      inputFiles: [],
      outputFiles: [],
      selectedFileContent: null,
      filePreviewLoading: false,
      filePreviewError: "",
      filePreviewRequestId: 0,
      plotImageUrl: "",
      plotLoading: false,
      plotError: "",
      plotOutput: "",
      plotUserGuidance: "",
      workflowContinuation: {
        eligible: false,
        reason: "",
        message: "",
        scientificVerification: null,
        nextStagePreview: null,
      },
      workflowContinuationLoading: false,
      plotForm: {
        file: "",
        x_axis: "0",
        y_axis: "1",
        flags: "-linY",
      },
    };
  },
  computed: {
    runId() {
      return this.$route.params.runId;
    },
    visualizationMode() {
      return this.$route.query.visualize === "1";
    },
    workflowSubtitle() {
      if (this.visualizationMode) {
        return "Workflow/Visualization";
      }
      if (this.presentation.subtitle === "Workflow/STGF") {
        return "Workflow/STGF";
      }
      return this.presentation.subtitle || "Workflow/Bound";
    },
    activeWorkflowStage() {
      return this.presentation.stages.find(stage => stage.state === "active")
          || this.presentation.stages.find(stage => stage.child_run_id === this.presentation.active_child_run_id)
          || null;
    },
    fileCatalogRunId() {
      if (
        this.presentation.mode === "workflow"
        && !this.isWorkflowChildRun
        && this.activeWorkflowStage
        && this.activeWorkflowStage.child_run_id
      ) {
        return this.activeWorkflowStage.child_run_id;
      }

      return this.runId;
    },
    workflowStatusLabel() {
      const state = this.presentation.workflow_state || "not_started";
      const labels = {
        not_started: "Not started",
        running: "Running",
        waiting: "Waiting for next step",
        complete: "Complete",
      };
      return labels[state] || state;
    },
    continuationPreviewItems() {
      const preview = this.workflowContinuation.nextStagePreview;
      if (!preview) {
        return [];
      }
      const stageLabels = {
        Data_Gen: "Data Generation",
        ePolyScat_Run: "ePolyScat Run",
        Analysis: "Visualization & Analysis",
      };

      return [
        {
          label: "Continue to",
          value: stageLabels[preview.nextStage]
              || this.formatStatusLabel(preview.nextStage),
          detail: preview.targetApplication || "",
        },
        {
          label: "Inherited input",
          value: this.continuationPreviewOutputName,
          detail: preview.inputFileName || "",
        },
      ];
    },
    continuationPreviewOutputName() {
      const preview = this.workflowContinuation.nextStagePreview;
      const outputFile = preview ? preview.outputFile : null;
      return outputFile
          ? this.fileDisplayName(outputFile)
          : "Manual selection required";
    },
    isWorkflowChildRun() {
      return Boolean(this.run && this.run.parentRunId);
    },
    workflowPanelItems() {
      if (this.visualizationMode) {
        const sourceStage = this.presentation.stages.find(stage => stage.id === "Analysis");
        return [
          { label: "Current stage", value: "Visualization" },
          {
            label: "Source step",
            value: sourceStage ? sourceStage.label : "Post-processing",
            detail: sourceStage ? sourceStage.application || "" : "",
          },
        ];
      }
      if (this.isWorkflowChildRun) {
        const items = [
          { label: "Step status", value: this.stepStatusLabel },
          { label: "Workflow position", value: this.workflowPositionLabel },
        ];

        if (this.nextWorkflowStage) {
          items.push({
            label: "Next step",
            value: this.nextWorkflowStage.label,
            detail: this.nextWorkflowStage.application || "",
          });
        }

        return items;
      }

      const items = [{ label: "Workflow status", value: this.workflowStatusLabel }];
      if (this.activeWorkflowStage) {
        items.push({
          label: "Current step",
          value: this.activeWorkflowStage.label,
          detail: this.activeWorkflowStage.application || "",
        });
      }

      return items;
    },
    stepStatusLabel() {
      if (!this.run) {
        return "";
      }

      return this.formatStatusLabel(
          this.run.displayStatus || this.run.status || this.run.workflowStepStatus
      );
    },
    workflowPositionLabel() {
      const stepOrder = this.workflowStepOrder();
      const totalSteps = this.presentation.stages.length || stepOrder;

      return stepOrder && totalSteps ? `Step ${stepOrder} of ${totalSteps}` : "Step";
    },
    nextWorkflowStage() {
      if (!this.isWorkflowChildRun) {
        return null;
      }

      const stepOrder = this.workflowStepOrder();
      return stepOrder ? this.presentation.stages[stepOrder] || null : null;
    },
    statusBadges() {
      if (!this.run) {
        return [];
      }

      const runStatus = this.run.displayStatus || this.run.status;
      const jobStatus = this.run.jobStatus;
      const badges = [];

      if (runStatus) {
        badges.push({ label: this.primaryStatusLabel, value: runStatus });
      }
      if (jobStatus && jobStatus !== runStatus) {
        badges.push({ label: "Job", value: jobStatus });
      }
      const verification = this.workflowContinuation.scientificVerification;
      if (verification && verification.status) {
        badges.push({
          label: "Science",
          value: this.formatStatusLabel(verification.status),
          title: verification.message || "",
        });
      }

      return badges;
    },
    primaryStatusLabel() {
      if (this.presentation.mode === "workflow") {
        if (this.isWorkflowChildRun) {
          return "Step";
        }
        return "Workflow";
      }
      return "Run";
    },
    runInputFiles() {
      if (!this.run || !this.run.inputs) {
        return [];
      }

      return this.run.inputs.reduce(
          (files, input) => files.concat(input.files || []),
          []
      );
    },
    showTargetStates() {
      return this.presentation.mode !== "workflow" || this.selectedFile !== "cross_sections";
    },
    targetInputColumns() {
      const actualInputFiles = this.resolvedInputFileNames;
      const selectedInputFiles = this.uniqueFileNames(this.presentationGroupFiles("Inputs"));

      if (actualInputFiles.length) {
        return this.chunkFileNames(actualInputFiles);
      }
      if (selectedInputFiles.length) {
        return this.chunkFileNames(selectedInputFiles);
      }

      return this.presentation.target_states.input_columns || [];
    },
    targetOutputFiles() {
      return this.resolvedOutputFiles;
    },
    activeRunApplicationForPlot() {
      if (this.activeWorkflowStage && this.activeWorkflowStage.application) {
        return this.activeWorkflowStage.application;
      }
      if (!this.run) {
        return "";
      }

      return this.run.utilityApplication
          || this.run.moduleApplication
          || this.run.workflowApplication
          || "";
    },
    isPlotCapableRun() {
      const stage = this.presentation.active_stage_id || (this.run ? this.run.workflowStage : "");
      const nonPlotApplications = ["Gaussian16", "OpenMolcas", "MoldenMerge"];

      if (stage === "Data_Gen" || stage === "data-generation" || stage === "Data_Generation") {
        return false;
      }

      return nonPlotApplications.indexOf(this.activeRunApplicationForPlot) < 0;
    },
    plottableOutputFiles() {
      return this.outputFiles.filter(file => (
        this.fileDataProductURI(file) && this.isPlottableFile(file)
      ));
    },
    showPlotPanel() {
      return this.visualizationMode
          || (this.isPlotCapableRun && this.hasPlottableOutputFiles);
    },
    resolvedInputFileNames() {
      return this.uniqueFileNames([
        ...this.inputFiles.map(file => this.fileDisplayName(file)),
        ...this.runInputFiles.map(file => this.fileDisplayName(file)),
        ...this.presentationGroupFiles("Inputs"),
      ]);
    },
    resolvedOutputFiles() {
      const inputNames = this.resolvedInputFileNames;
      const outputFileNames = this.outputFiles.map(file => this.fileDisplayName(file));
      const viewableNames = this.viewables.map(file => this.fileDisplayName(file));

      return this.uniqueFileNames([
        ...outputFileNames.filter(filename => inputNames.indexOf(filename) < 0),
        ...viewableNames.filter(filename => inputNames.indexOf(filename) < 0),
      ]);
    },
    resolvedFileGroups() {
      const inputNames = this.resolvedInputFileNames;
      const outputNames = this.resolvedOutputFiles;
      const groups = [];

      if (inputNames.length) {
        groups.push({ label: "Inputs", files: inputNames });
      }
      if (outputNames.length) {
        groups.push({ label: "Outputs", files: outputNames });
      }

      return groups.length ? groups : this.presentation.file_groups;
    },
    plotFileOptions() {
      const files = this.plottableOutputFiles;

      if (!files.length) {
        return [{ value: "", text: "No plottable output files yet", disabled: true }];
      }

      return files.map(file => {
        const contract = this.plotContractForFile(file);
        const dimension = contract && contract.dimension
          ? ` (${contract.dimension}D)`
          : "";

        return {
          value: this.fileDisplayName(file),
          text: `${this.fileDisplayName(file)}${dimension}`,
        };
      });
    },
    hasPlottableOutputFiles() {
      return this.plottableOutputFiles.length > 0;
    },
    canCreatePlot() {
      const outputFile = this.plottableOutputFileForName(this.plotForm.file);

      return Boolean(
          outputFile
          && this.fileDataProductURI(outputFile)
          && this.plotForm.x_axis
          && this.plotForm.y_axis
      );
    },
    codePreview() {
      if (this.filePreviewLoading) {
        return `# ${this.selectedFile}\nLoading file preview...`;
      }
      if (this.filePreviewError) {
        return `# ${this.selectedFile}\n${this.filePreviewError}`;
      }
      if (this.selectedFileContent !== null) {
        return this.selectedFileContent;
      }
      if (this.presentation.code && this.selectedFile === this.presentation.selected_file) {
        return this.presentation.code;
      }

      return [
        `# ${this.selectedFile}`,
        "No preview available for this file yet.",
      ].join("\n");
    },
  },
  watch: {
    selectedFile() {
      this.fetchSelectedFileContent();
    },
    runId() {
      this.refreshRun();
    },
    visualizationMode(enabled) {
      if (!enabled) {
        return;
      }
      this.initializePlotForm();
      if (this.plotForm.file) {
        this.selectedFile = this.plotForm.file;
      }
    },
    "plotForm.file"(filename) {
      const outputFile = this.plottableOutputFileForName(filename);
      this.applyPlotContractForFile(outputFile);
      if (filename && this.visualizationMode) {
        this.selectedFile = filename;
      }
    },
  },
  methods: {
    async refreshRun() {
      try {
        this.run = await RunService.fetchRun({ runId: this.runId });
        this.presentation = this.normalizePresentation(this.run.presentation);
        await this.loadRunFileCatalog();
        await this.loadWorkflowContinuation();
        this.initializePlotForm();

        const nextSelectedFile = (this.visualizationMode && this.plotForm.file)
            || this.presentation.selected_file
            || this.firstResolvedFile()
            || this.firstPresentationFile();
        if (this.selectedFile === nextSelectedFile) {
          await this.fetchSelectedFileContent();
        } else {
          this.selectedFile = nextSelectedFile;
        }
      } catch (error) {
        eventBus.$emit("error", { name: `Error while trying to fetch run with id: ${this.runId}`, error });
      }
    },
    normalizePresentation(presentation) {
      const normalized = {
        mode: "module",
        subtitle: "Modules/EPOLYSCAT_DMAT",
        file_groups: [],
        target_states: { input_columns: [], output_files: [] },
        parameters: [],
        plot: { file: "", x_axis: "0", y_axis: "1", flags: "-linY" },
        plottable_file_names: [],
        stages: [],
        code: "",
        ...presentation,
      };

      normalized.file_groups = normalized.file_groups || [];
      normalized.target_states = {
        input_columns: [],
        output_files: [],
        ...(normalized.target_states || {}),
      };
      normalized.plot = {
        file: "",
        x_axis: "0",
        y_axis: "1",
        flags: "-linY",
        ...(normalized.plot || {}),
      };
      normalized.stages = normalized.stages || [];
      normalized.plottable_file_names = normalized.plottable_file_names || [];

      return normalized;
    },
    async loadWorkflowContinuation() {
      this.workflowContinuation = {
        eligible: false,
        reason: "",
        message: "",
        scientificVerification: null,
        nextStagePreview: null,
      };
      try {
        this.workflowContinuation = await RunService.fetchWorkflowContinuation({
          runId: this.runId,
        });
      } catch (error) {
        this.workflowContinuation = {
          eligible: false,
          reason: "status_unavailable",
          message: "Workflow continuation is unavailable.",
          scientificVerification: null,
          nextStagePreview: null,
        };
      }
    },
    async continueInWorkflow() {
      if (!this.workflowContinuation.eligible || this.workflowContinuationLoading) {
        return;
      }

      this.workflowContinuationLoading = true;
      try {
        const continuation = await RunService.continueWorkflow({ runId: this.runId });
        await this.$router.push(
            `/runs/new?workflowChildRunId=${continuation.nextChildRunId}`
            + `&workflowParentRunId=${continuation.workflowParentRunId}`
            + `&withOutputsFrom=${continuation.sourceRunId}`
        );
      } catch (error) {
        eventBus.$emit("error", {
          name: "Error while trying to continue this run in a workflow",
          error,
        });
        await this.loadWorkflowContinuation();
      } finally {
        this.workflowContinuationLoading = false;
      }
    },
    async loadRunFileCatalog() {
      const [inputFiles, outputFiles] = await Promise.all([
        PlotService.getInputFiles({ runId: this.fileCatalogRunId }).catch(() => []),
        InputService.fetchOutputs(this.fileCatalogRunId).catch(() => []),
      ]);

      this.inputFiles = inputFiles || [];
      this.outputFiles = outputFiles || [];
      this.viewables = this.outputFiles.filter(file => file.viewable);
    },
    initializePlotForm() {
      const currentOutput = this.plottableOutputFileForName(this.plotForm.file || this.presentation.plot.file);
      const firstOutput = this.plottableOutputFiles[0];
      const selectedOutput = currentOutput && this.fileDataProductURI(currentOutput)
          ? currentOutput
          : firstOutput || null;
      const plotFile = selectedOutput ? this.fileDisplayName(selectedOutput) : "";

      this.plotForm = {
        file: plotFile,
        x_axis: this.presentation.plot.x_axis || "0",
        y_axis: this.presentation.plot.y_axis || "1",
        flags: this.presentation.plot.flags || "",
      };
      this.applyPlotContractForFile(selectedOutput);
      this.plotImageUrl = "";
      this.plotError = "";
      this.plotOutput = "";
      this.plotUserGuidance = "";
    },
    async fetchSelectedFileContent() {
      const requestId = this.filePreviewRequestId + 1;
      this.filePreviewRequestId = requestId;
      this.selectedFileContent = null;
      this.filePreviewError = "";

      if (!this.selectedFile) {
        return;
      }

      this.filePreviewLoading = true;
      try {
        const inputFile = this.inputFileForName(this.selectedFile);
        const outputFile = this.outputFileForName(this.selectedFile);
        if (outputFile && outputFile.viewable === false) {
          this.filePreviewError = "This file is binary and cannot be previewed as text.";
          return;
        }
        const content = inputFile && this.fileDataProductURI(inputFile)
            ? await InputService.fetchFileContents(inputFile)
            : outputFile && this.fileDataProductURI(outputFile)
            ? await InputService.fetchFileContents(outputFile)
            : await RunService.fetchViewableContent({
              runId: this.fileCatalogRunId,
              filename: this.selectedFile,
            });

        if (requestId !== this.filePreviewRequestId) {
          return;
        }

        if (content === null || content === undefined || content === "") {
          this.filePreviewError = "This file is not a text file or is not available yet.";
        } else {
          this.selectedFileContent = content;
        }
      } catch (error) {
        if (requestId === this.filePreviewRequestId) {
          this.filePreviewError = "No preview is available for this file yet.";
        }
      } finally {
        if (requestId === this.filePreviewRequestId) {
          this.filePreviewLoading = false;
        }
      }
    },
    async createPlot() {
      const outputFile = this.plottableOutputFileForName(this.plotForm.file);

      if (!outputFile || !this.fileDataProductURI(outputFile)) {
        this.plotError = "No plottable output file is available yet.";
        return;
      }

      this.plotLoading = true;
      this.plotError = "";
      this.plotOutput = "";
      this.plotUserGuidance = "";
      try {
        const result = await PlotService.plotSelectedRuns({
          runIds: [this.fileCatalogRunId],
          plotfiles: [{
            dataProductURI: this.fileDataProductURI(outputFile),
            prefix: "",
          }],
          xAxis: this.plotForm.x_axis,
          yAxis: this.plotForm.y_axis,
          flags: this.plotForm.flags,
        });

        this.plotImageUrl = result.plotImageUrl || "";
        this.plotOutput = result.output || "";
        this.plotUserGuidance = result.userGuidance || "";
        if (!this.plotImageUrl && !this.plotOutput && !this.plotUserGuidance) {
          this.plotError = "The plot command finished without a plot image.";
        }
      } catch (error) {
        this.plotImageUrl = "";
        this.plotError = "Plot failed for the selected output file.";
        eventBus.$emit("error", { name: "Error while trying to plot this run", error });
      } finally {
        this.plotLoading = false;
      }
    },
    firstPresentationFile() {
      const firstGroup = this.presentation.file_groups[0];
      if (firstGroup && firstGroup.files.length) {
        return firstGroup.files[0];
      }
      return "target";
    },
    firstResolvedFile() {
      const firstGroup = this.resolvedFileGroups[0];

      if (firstGroup && firstGroup.files.length) {
        return firstGroup.files[0];
      }

      return "";
    },
    presentationGroupFiles(label) {
      const group = this.presentation.file_groups.find(fileGroup => fileGroup.label === label);

      return group ? group.files : [];
    },
    uniqueFileNames(filenames) {
      return [...new Set(filenames.filter(filename => Boolean(filename)))];
    },
    chunkFileNames(filenames) {
      const filteredNames = this.uniqueFileNames(filenames);
      const chunkSize = Math.max(1, Math.ceil(filteredNames.length / 3));
      const chunks = [];

      for (let index = 0; index < filteredNames.length; index += chunkSize) {
        chunks.push(filteredNames.slice(index, index + chunkSize));
      }

      return chunks;
    },
    fileDisplayName(file) {
      if (typeof file === "string") {
        return file;
      }

      return file.filename || file.name || file.fileName || "";
    },
    fileDataProductURI(file) {
      if (!file) {
        return "";
      }

      return file.dataProductURI || file["data-product-uri"] || file.data_product_uri || "";
    },
    outputFileForName(filename) {
      return this.outputFiles.find(file => this.fileDisplayName(file) === filename);
    },
    plottableOutputFileForName(filename) {
      return this.plottableOutputFiles.find(file => this.fileDisplayName(file) === filename);
    },
    isPlottableFile(file) {
      return Boolean(file && (file.plottable === true || file.plot_contract));
    },
    plotContractForFile(file) {
      if (!file) {
        return null;
      }

      return file.plot_contract || file.plotContract || null;
    },
    applyPlotContractForFile(file) {
      const contract = this.plotContractForFile(file);
      if (!contract) {
        return;
      }

      this.plotForm.x_axis = contract.x_axis || "0";
      this.plotForm.y_axis = contract.y_axes || contract.y_axis || "1";
      this.plotForm.flags = contract.flags || "";
    },
    inputFileForName(filename) {
      return this.runInputFiles.find(file => this.fileDisplayName(file) === filename);
    },
    workflowStepOrder() {
      const runId = this.run ? this.run.id : null;
      const presentationIndex = this.presentation.stages.findIndex(stage => (
        stage.child_run_id === runId
        || (stage.child_run_ids || []).indexOf(runId) >= 0
      ));
      if (presentationIndex >= 0) {
        return presentationIndex + 1;
      }
      if (this.run && this.run.workflowStepOrder) {
        return this.run.workflowStepOrder;
      }

      const activeIndex = this.presentation.stages.findIndex(stage => stage === this.activeWorkflowStage);
      return activeIndex >= 0 ? activeIndex + 1 : 0;
    },
    workflowStageIsConfigurable(stage, index) {
      if (!stage || stage.local_only || !stage.child_run_id || !this.run || this.isWorkflowChildRun) {
        return false;
      }
      if ((stage.status || "pending") !== "pending") {
        return false;
      }

      return this.presentation.stages.slice(0, index).every(previousStage =>
        previousStage.state === "complete"
        || previousStage.status === "complete"
        || previousStage.status === "imported"
        || previousStage.status === "not_included"
      );
    },
    workflowStageRoute(stage, index) {
      if (!stage) {
        return null;
      }
      if (stage.local_only) {
        return stage.status === "ready" && stage.source_child_run_id
          ? { path: `/runs/${stage.source_child_run_id}`, query: { visualize: "1" } }
          : null;
      }
      if (!stage.child_run_id) {
        return null;
      }
      if (this.workflowStageIsConfigurable(stage, index)) {
        return this.workflowConfigureStepLink(stage, index);
      }
      return `/runs/${stage.child_run_id}`;
    },
    workflowConfigureStepLink(stage, index) {
      const query = [
        `workflowChildRunId=${stage.child_run_id}`,
        `workflowParentRunId=${this.run.id}`,
      ];
      const previousStage = this.presentation.stages
          .slice(0, index)
          .reverse()
          .find(item => item.child_run_id);

      if (previousStage) {
        query.push(`withOutputsFrom=${previousStage.child_run_id}`);
      }

      return `/runs/new?${query.join("&")}`;
    },
    workflowStageStatusText(stage) {
      if (!stage) {
        return "";
      }
      if (stage.status === "imported") {
        return "Imported";
      }
      if (stage.status === "not_included") {
        return "Not included";
      }
      if (stage.status === "ready") {
        return "Ready";
      }
      return "";
    },
    formatStatusLabel(status) {
      if (!status) {
        return "";
      }

      return status
          .toString()
          .toLowerCase()
          .replace(/_/g, " ")
          .replace(/\b\w/g, character => character.toUpperCase());
    },
  },
  mounted() {
    this.refreshRun();
  },
};
</script>

<style scoped>
.view-run-page {
  background: #ffffff;
  color: #000000;
  min-height: 100%;
}

.view-run-content {
  max-width: 1200px;
  padding: 18px 32px 90px 44px;
}

.view-run-breadcrumb {
  font-size: 16px;
  margin-bottom: 18px;
}

.view-run-heading {
  align-items: flex-start;
  display: flex;
  justify-content: space-between;
  margin-bottom: 24px;
}

.view-run-heading h1 {
  font-size: 28px;
  font-weight: 700;
  line-height: 1.2;
  margin: 0 0 6px;
}

.view-run-subtitle {
  color: #555555;
  font-size: 17px;
  font-weight: 700;
}

.view-run-status {
  display: flex;
  gap: 8px;
}

.view-run-heading-actions {
  align-items: flex-end;
  display: flex;
  flex-direction: column;
  gap: 10px;
}

.view-run-status span {
  background: #f1f5f7;
  border-radius: 6px;
  color: #226597;
  font-size: 13px;
  font-weight: 700;
  padding: 6px 10px;
}

.workflow-continuation-preview {
  border-bottom: 1px solid #d8dde5;
  border-top: 1px solid #d8dde5;
  display: grid;
  gap: 24px;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  margin: -4px 0 28px;
  max-width: 1017px;
  padding: 14px 0;
}

.workflow-continuation-preview > div {
  display: flex;
  flex-direction: column;
  gap: 2px;
  min-width: 0;
}

.workflow-continuation-preview strong,
.workflow-continuation-preview .workflow-status-application {
  overflow-wrap: anywhere;
}

.workflow-stepper {
  align-items: center;
  display: flex;
  justify-content: center;
  margin-bottom: 30px;
}

.workflow-step {
  align-items: center;
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

.workflow-step-label {
  color: #222222;
}

.workflow-step-link {
  text-decoration: none;
}

.workflow-step-link:hover {
  color: #226597;
  text-decoration: underline;
}

.workflow-step-configure-link {
  color: #226597;
}

.workflow-step-application {
  color: #666666;
  font-size: 12px;
  font-weight: 500;
}

.workflow-step-copy {
  display: flex;
  flex-direction: column;
  gap: 1px;
}

.workflow-step-state {
  color: #666666;
  font-size: 12px;
  font-weight: 600;
}

.workflow-status-panel {
  align-items: center;
  background: #f7f8fa;
  border: 1px solid #d8dde5;
  border-radius: 6px;
  display: flex;
  gap: 18px;
  justify-content: space-between;
  margin: -12px 0 30px;
  max-width: 1017px;
  padding: 14px 16px;
}

.workflow-status-panel > div {
  display: flex;
  flex-direction: column;
  gap: 2px;
}

.workflow-status-label {
  color: #666666;
  font-size: 12px;
  font-weight: 700;
  text-transform: uppercase;
}

.workflow-status-application {
  color: #555555;
  font-size: 13px;
  font-weight: 600;
}

.view-run-layout {
  display: grid;
}

.view-run-layout.with-plot-panel {
  display: grid;
  gap: 28px;
  grid-template-columns: minmax(0, 1fr) 360px;
}

.view-run-layout:not(.with-plot-panel) .view-run-main {
  justify-self: center;
  max-width: 760px;
  width: 100%;
}

.file-control-row {
  align-items: center;
  display: flex;
  gap: 12px;
  margin-bottom: 0;
}

.file-control-row label {
  font-weight: 700;
  margin: 0;
}

.file-selector {
  border: 1px solid #d6d6d6;
  border-radius: 6px;
  height: 38px;
  max-width: 320px;
  min-width: 220px;
  padding: 0 10px;
}

.target-states-matrix {
  align-items: start;
  display: grid;
  gap: 22px;
  grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
  margin-bottom: 22px;
}

.target-column h2,
.view-run-resource h2 {
  font-size: 17px;
  font-weight: 700;
  margin: 0 0 12px;
}

.target-input-grid {
  display: grid;
  gap: 10px;
  grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
}

.target-input-column,
.target-output-list {
  align-items: start;
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
}

.target-empty-state {
  background: #f7f8fa;
  border: 1px dashed #cfd6df;
  border-radius: 6px;
  color: #666666;
  font-size: 13px;
  font-weight: 700;
  padding: 12px;
}

.file-chip,
.module-tab {
  background: #eeeeee;
  border: 0;
  border-radius: 6px;
  color: #000000;
  font-size: 15px;
  font-weight: 700;
  min-height: 36px;
  min-width: 0;
  overflow-wrap: anywhere;
  padding: 8px 12px;
}

.file-chip.active,
.module-tab.active {
  background: #226597;
  color: #ffffff;
}

.module-parameters {
  border: 1px solid #dddddd;
  border-radius: 6px;
  margin-bottom: 22px;
  max-width: 760px;
}

.module-tabs {
  border-bottom: 1px solid #dddddd;
  display: flex;
  gap: 8px;
  padding: 14px;
}

.module-parameter-grid {
  display: grid;
  gap: 10px;
  padding: 16px;
}

.module-parameter-row {
  display: grid;
  gap: 12px;
  grid-template-columns: 220px 120px;
}

.run-code-viewer {
  border: 1px solid #dddddd;
  border-radius: 6px;
  margin-bottom: 28px;
  max-width: 760px;
}

.code-viewer-toolbar {
  align-items: center;
  border-bottom: 1px solid #dddddd;
  display: flex;
  gap: 14px;
  font-weight: 700;
  justify-content: space-between;
  padding: 10px 14px;
}

.code-viewer-toolbar > span {
  min-width: 0;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}

.code-viewer-toolbar .file-control-row {
  flex: 0 0 auto;
}

.run-code-viewer pre {
  background: #fafafa;
  color: #1d1d1d;
  font-size: 13px;
  line-height: 1.55;
  margin: 0;
  min-height: 180px;
  overflow: auto;
  padding: 16px;
  white-space: pre-wrap;
}

.view-run-resource {
  max-width: 760px;
}

.plot-panel {
  border-left: 1px solid #dddddd;
  display: grid;
  gap: 14px;
  padding-left: 24px;
}

.plot-panel label {
  display: grid;
  font-size: 13px;
  font-weight: 700;
  gap: 6px;
}

.plot-canvas {
  align-items: center;
  background: #ffffff;
  border: 1px solid #dddddd;
  border-radius: 6px;
  display: flex;
  height: 220px;
  justify-content: center;
  overflow: hidden;
}

.plot-canvas img {
  height: 100%;
  object-fit: contain;
  width: 100%;
}

.plot-placeholder {
  color: #666666;
  font-size: 13px;
  font-weight: 700;
  line-height: 1.4;
  padding: 18px;
  text-align: center;
}

.plot-output {
  background: #f7f8fa;
  border: 1px solid #d8dde5;
  border-radius: 6px;
  color: #1d1d1d;
  font-size: 12px;
  line-height: 1.45;
  margin: 0;
  max-height: 160px;
  overflow: auto;
  padding: 10px;
  white-space: pre-wrap;
}

.plot-guidance {
  color: #666666;
  font-size: 13px;
  font-weight: 600;
  line-height: 1.4;
}

@media (max-width: 980px) {
  .view-run-content {
    padding: 16px 18px 70px;
  }

  .view-run-heading,
  .view-run-layout.workflow-layout {
    display: block;
  }

  .plot-panel {
    border-left: 0;
    border-top: 1px solid #dddddd;
    margin-top: 28px;
    padding-left: 0;
    padding-top: 24px;
  }

  .target-states-matrix,
  .target-input-grid,
  .workflow-continuation-preview {
    grid-template-columns: 1fr;
  }
}
</style>
