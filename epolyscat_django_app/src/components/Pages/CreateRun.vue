<template>
  <div class="new-run-page">
    <main class="new-run-content" ref="newRunContent">
      <span
          id="new-run-top-anchor"
          ref="newRunTopAnchor"
          class="new-run-top-anchor"
          tabindex="-1"
          aria-hidden="true"
      />
      <b-breadcrumb class="new-run-breadcrumb">
        <template v-if="experimentIdFromQueryString">
          <b-breadcrumb-item to="/experiments">Experiments</b-breadcrumb-item>
          <b-breadcrumb-item v-if="experiment" :to="`/runs?experimentId=${experiment.experimentId}`">
            {{ experiment.name }}
          </b-breadcrumb-item>
        </template>
        <template v-else-if="viewId">
          <b-breadcrumb-item to="/views">Views</b-breadcrumb-item>
          <b-breadcrumb-item v-if="view" :to="`/runs?viewId=${view.viewId}`">
            {{ view.name }}
          </b-breadcrumb-item>
        </template>
        <b-breadcrumb-item v-else to="/runs">Experiments</b-breadcrumb-item>

        <template v-if="cloneRunId">
          <b-breadcrumb-item :to="runLink" v-if="cloneRunOriginalName">
            {{ cloneRunOriginalName }}
          </b-breadcrumb-item>
          <b-breadcrumb-item :to="newRunLink">Clone</b-breadcrumb-item>
        </template>
        <b-breadcrumb-item v-else :to="newRunLink">New Run</b-breadcrumb-item>
      </b-breadcrumb>

      <header class="new-run-heading">
        <h1 class="new-run-title">{{ pageTitle }}</h1>
      </header>

      <section class="run-type-selector-section">
        <template v-if="selectedRunType !== 'workflow'">
          <div class="run-type-selection-grid">
            <div class="run-selection-column run-type-column">
              <h2>Run Type</h2>
              <button
                  v-for="runType in runTypeOptions"
                  :key="runType.id"
                  type="button"
                  class="run-selection-option"
                  :class="{ active: selectedRunType === runType.id }"
                  v-on:click="selectRunType(runType.id)"
              >
                <span class="run-selection-label">{{ runType.label }}</span>
                <span class="run-selection-description">{{ runType.description }}</span>
              </button>
            </div>

            <div class="run-selection-column application-column">
              <h2>{{ activeRunTypeTitle }}</h2>
              <button
                  v-for="application in activeRunTypeApplications"
                  :key="application.id"
                  type="button"
                  class="run-selection-option"
                  :class="{ active: selectedApplicationId === application.id }"
                  v-on:click="selectRunApplication(application.id)"
              >
                <span class="run-selection-label">{{ application.label }}</span>
                <span class="run-selection-description">{{ application.description }}</span>
              </button>
            </div>

            <div class="run-selection-column required-files-column">
              <h2>Required Files</h2>
              <div class="required-file-list">
                <div
                    v-for="file in activeRequiredFiles"
                    :key="file.name"
                    class="required-file-row"
                    role="button"
                    tabindex="0"
                    v-on:click="selectRequiredFile(file.name)"
                    v-on:keydown.enter.prevent="selectRequiredFile(file.name)"
                >
                  <span class="required-file-name">{{ file.label }}</span>
                  <span class="required-file-description">{{ file.description }}</span>
                  <span
                      class="required-file-status"
                      :class="uploadedRequiredFiles(file.name).length > 0 ? 'required-file-status-ready' : 'required-file-status-missing'"
                  >
                    {{ uploadedRequiredFiles(file.name).length > 0 ? "Ready" : "Missing" }}
                  </span>
                  <span
                      v-if="uploadedRequiredFiles(file.name).length > 0"
                      class="required-file-selected-list"
                  >
                    {{ uploadedRequiredFiles(file.name).map(uploadedFile => uploadedFile.name).join(", ") }}
                  </span>
                </div>
                <div v-if="activeRequiredFiles.length === 0" class="required-file-empty">
                  No required files
                </div>
              </div>
            </div>
          </div>
        </template>

        <template v-else>
          <div class="run-type-switcher-row" role="tablist" aria-label="Run Type">
            <button
                v-for="runType in runTypeOptions"
                :key="runType.id"
                type="button"
                class="run-type-switch-option"
                :class="{ active: selectedRunType === runType.id }"
                :aria-selected="selectedRunType === runType.id"
                v-on:click="selectRunType(runType.id)"
            >
              <span class="run-type-switch-label">{{ runType.label }}</span>
              <span class="run-type-switch-description">{{ runType.description }}</span>
            </button>
          </div>

          <nav class="workflow-stepper" aria-label="Run workflow">
            <button
                v-for="(application, index) in activeRunTypeApplications"
                :key="application.id"
                type="button"
                class="workflow-step"
                :class="{ active: selectedApplicationId === application.id, disabled: application.localOnly }"
                :disabled="application.localOnly"
                :aria-current="selectedApplicationId === application.id ? 'step' : null"
                v-on:click="selectRunApplication(application.id)"
            >
              <span class="workflow-step-index">{{ index + 1 }}</span>
              <span class="workflow-step-label">{{ application.label }}</span>
            </button>
          </nav>

          <div class="workflow-selection-grid">
            <div class="run-selection-column application-column">
              <h2>{{ activeWorkflowApplicationTitle }}</h2>
              <div
                  v-if="activeWorkflowApplications.length > 0"
                  class="workflow-stage-application-panel"
              >
                <div class="workflow-stage-application-options">
                  <button
                      v-for="application in activeWorkflowApplications"
                      :key="application.id"
                      type="button"
                      class="workflow-stage-application-option"
                      :class="{ active: activeWorkflowApplicationId === application.id }"
                      v-on:click="selectWorkflowApplication(application.id)"
                  >
                    <span class="workflow-stage-application-label">{{ application.label }}</span>
                    <span class="workflow-stage-application-description">{{ application.description }}</span>
                  </button>
                </div>
              </div>
            </div>

            <div class="run-selection-column required-files-column">
              <h2>Required Files</h2>
              <div class="required-file-list">
                <div
                    v-for="file in activeRequiredFiles"
                    :key="file.name"
                    class="required-file-row"
                    role="button"
                    tabindex="0"
                    v-on:click="selectRequiredFile(file.name)"
                    v-on:keydown.enter.prevent="selectRequiredFile(file.name)"
                >
                  <span class="required-file-name">{{ file.label }}</span>
                  <span class="required-file-description">{{ file.description }}</span>
                  <span
                      class="required-file-status"
                      :class="uploadedRequiredFiles(file.name).length > 0 ? 'required-file-status-ready' : 'required-file-status-missing'"
                  >
                    {{ uploadedRequiredFiles(file.name).length > 0 ? "Ready" : "Missing" }}
                  </span>
                  <span
                      v-if="uploadedRequiredFiles(file.name).length > 0"
                      class="required-file-selected-list"
                  >
                    {{ uploadedRequiredFiles(file.name).map(uploadedFile => uploadedFile.name).join(", ") }}
                  </span>
                </div>
                <div v-if="activeRequiredFiles.length === 0" class="required-file-empty">
                  No required files
                </div>
              </div>
            </div>
          </div>
        </template>
      </section>

      <section class="new-run-descriptions">
        <div class="new-run-description-header">
          <h2>Descriptions</h2>
          <b-button
              variant="link"
              class="new-run-description-toggle"
              v-on:click="showDescriptions = !showDescriptions"
          >
            <span v-if="!showDescriptions">(Show)</span>
            <span v-else>(Hide)</span>
          </b-button>
        </div>
        <b-collapse
            id="new-run-description-list"
            class="new-run-description-list"
            v-model="showDescriptions"
        >
          <div
              v-for="description in newRunDescriptions"
              :key="description.name"
              class="new-run-description-row"
          >
            <b>{{ description.name }}:</b> {{ description.text }}
          </div>
        </b-collapse>
      </section>

      <b-alert
          v-if="workflowOutputBindingMessage"
          show
          variant="warning"
          class="workflow-output-binding-warning"
      >
        {{ workflowOutputBindingMessage }}
      </b-alert>

      <section v-if="hasRequiredFiles" class="input-files-panel" ref="inputFilesPanel">
        <h2>Input Files</h2>
        <div class="input-file-tabs" role="tablist" aria-label="Input files">
          <button
              v-for="file in inputFiles"
              :key="file.id"
              type="button"
              class="input-file-tab"
              :class="{ active: selectedInputFile === file.id }"
              :aria-selected="selectedInputFile === file.id"
              v-on:click="selectInputFile(file.id)"
          >
            {{ file.label }}
          </button>
        </div>
        <div v-if="activeInputFile.inputFileName" class="input-file-actions">
          <b-button
              variant="dark"
              class="input-file-storage-button"
              v-b-modal.new-run-input-file-storage
          >
            Select file from storage
          </b-button>
          <span class="input-file-or">OR</span>
          <label class="input-file-computer-button" for="new-run-input-file-computer">
            Select file from computer
          </label>
          <input
              id="new-run-input-file-computer"
              class="input-file-computer-input"
              type="file"
              :multiple="activeInputAllowsMultiple"
              v-on:change="onInputComputerFilesSelected"
          >
        </div>
        <div v-if="activeInputFiles.length > 0" class="input-file-selected-list">
          <span
              v-for="file in activeInputFiles"
              :key="file.name"
              class="input-file-selected"
          >
            <span class="input-file-selected-name">{{ file.name }}</span>
            <button
                type="button"
                class="input-file-remove"
                :aria-label="`Remove ${file.name}`"
                v-on:click="removeInputFile(file)"
            >
              <b-icon icon="x" aria-hidden="true" />
            </button>
          </span>
        </div>
        <UserStorage
            id="new-run-input-file-storage"
            selectionType="files"
            :canSelectMultiple="activeInputAllowsMultiple"
            modalTitle="Select required input file"
            :onlyShowEditable="false"
            @filesSelected="onInputStorageFilesSelected"
        />
      </section>

      <section v-if="hasRequiredFiles && isEPolyScatScriptInput" class="data-entry-section">
        <div class="section-title-row">
          <h2>Data Entry</h2>
          <div class="view-mode-switcher" aria-label="Data entry view mode">
            <button
                type="button"
                class="view-mode-button"
                :class="{ active: dataEntryViewMode === 'file' }"
                v-on:click="selectDataEntryViewMode('file')"
                title="File view"
            >
              <b-icon icon="file-earmark-text" />
            </button>
            <button
                type="button"
                class="view-mode-button"
                :class="{ active: dataEntryViewMode === 'table' }"
                v-on:click="selectDataEntryViewMode('table')"
                title="Table view"
            >
              <b-icon icon="table" />
            </button>
          </div>
        </div>

        <div class="data-entry-card" :class="dataEntryCardClasses">
          <template v-if="dataEntryViewMode === 'table'">
            <div class="data-entry-tabs" role="tablist" aria-label="Data entry tabs">
              <button
                  v-for="section in dataEntrySections"
                  :key="section.id"
                  type="button"
                  class="data-entry-tab"
                  :class="{ active: selectedDataEntryTab === section.id }"
                  v-on:click="selectDataEntryTab(section.id)"
              >
                {{ section.label }}
              </button>
            </div>

            <div class="data-entry-body data-entry-body-table">
              <div class="data-entry-editor-grid">
                <section class="data-entry-editor-panel">
                  <div class="data-entry-section-heading">
                    <h3>{{ activeDataEntryTabLabel }}</h3>
                    <span>{{ activeDataEntrySection.recordGroup }}</span>
                  </div>

                  <div
                      v-if="activeDataEntrySection.type === 'fields'"
                      class="manual-field-rows"
                  >
                    <div
                        v-for="field in activeDataEntryFields"
                        :key="field.key"
                        class="manual-field-row"
                    >
                      <label :for="`data-entry-${field.key}`">{{ field.label }}</label>
                      <span class="manual-field-record">{{ field.record }}</span>
                      <b-form-select
                          v-if="field.inputType === 'select'"
                          :id="`data-entry-${field.key}`"
                          v-model="dataEntryValues[field.key]"
                          :options="fieldSelectOptions(field)"
                          class="manual-field-input"
                          v-on:input="onTableFieldInput(field.key, $event)"
                      />
                      <b-form-input
                          v-else
                          :id="`data-entry-${field.key}`"
                          v-model="dataEntryValues[field.key]"
                          :type="field.inputType === 'number' ? 'number' : 'text'"
                          :step="field.step || null"
                          :placeholder="dataEntryPlaceholder(field.key)"
                          class="manual-field-input"
                          v-on:input="onTableFieldInput(field.key, $event)"
                      />
                    </div>
                  </div>

                  <div
                      v-else-if="activeDataEntrySection.type === 'outputs'"
                      class="output-record-grid"
                  >
                    <div
                        v-for="output in outputDefinitions"
                        :key="output.fileType"
                        class="output-record-row"
                    >
                      <span class="output-record-type">{{ output.fileType }}</span>
                      <span class="output-record-description">{{ output.description }}</span>
                      <b-form-input
                          :id="`output-${output.valueKey}`"
                          v-model="dataEntryValues[output.valueKey]"
                          :placeholder="dataEntryPlaceholder(output.valueKey)"
                          class="manual-field-input"
                          v-on:input="onTableFieldInput(output.valueKey, $event)"
                      />
                      <span class="output-record-extension">{{ output.extension }}</span>
                    </div>
                  </div>

                  <div
                      v-else-if="activeDataEntrySection.type === 'sequence'"
                      class="ordered-sequence-editor"
                  >
                    <div class="ordered-sequence-toolbar">
                      <b-form-select
                          id="sequence-node-type"
                          v-model="sequenceNodeType"
                          :options="sequenceNodeTypeOptions"
                          aria-label="Record or command type"
                          class="ordered-sequence-select"
                          v-on:input="onSequenceNodeTypeInput"
                      />
                      <b-form-select
                          id="sequence-node-label"
                          v-model="sequenceNodeLabel"
                          :options="sequenceNodeLabelOptions"
                          aria-label="Record or command name"
                          class="ordered-sequence-select ordered-sequence-label-select"
                      />
                      <b-button
                          type="button"
                          variant="primary"
                          class="ordered-sequence-add"
                          v-on:click="appendDataEntrySequenceNode"
                      >
                        <b-icon icon="plus" aria-hidden="true" />
                        Add
                      </b-button>
                    </div>

                    <div class="ordered-sequence-list">
                      <div
                          v-for="(item, sequenceIndex) in dataEntrySequenceNodes"
                          :key="`${item.nodeIndex}-${item.type}-${item.label}-${item.occurrence}`"
                          class="ordered-sequence-row"
                      >
                        <div class="ordered-sequence-row-header">
                          <span class="ordered-sequence-index">{{ sequenceIndex + 1 }}</span>
                          <div class="ordered-sequence-identity">
                            <strong>{{ item.label }}</strong>
                            <span>{{ sequenceNodeTypeText(item.type) }} - occurrence {{ item.occurrence }}</span>
                          </div>
                          <div class="ordered-sequence-actions">
                            <button
                                type="button"
                                class="ordered-sequence-action"
                                :disabled="sequenceIndex === 0"
                                title="Move up"
                                aria-label="Move up"
                                v-on:click="moveDataEntrySequenceNode(item.nodeIndex, 'up')"
                            >
                              <b-icon icon="chevron-up" aria-hidden="true" />
                            </button>
                            <button
                                type="button"
                                class="ordered-sequence-action"
                                :disabled="sequenceIndex === dataEntrySequenceNodes.length - 1"
                                title="Move down"
                                aria-label="Move down"
                                v-on:click="moveDataEntrySequenceNode(item.nodeIndex, 'down')"
                            >
                              <b-icon icon="chevron-down" aria-hidden="true" />
                            </button>
                            <button
                                type="button"
                                class="ordered-sequence-action ordered-sequence-remove"
                                title="Remove"
                                aria-label="Remove"
                                v-on:click="removeDataEntrySequenceNode(item.nodeIndex)"
                            >
                              <b-icon icon="x" aria-hidden="true" />
                            </button>
                          </div>
                        </div>
                        <b-form-textarea
                            :value="item.raw"
                            :rows="item.multiline ? 5 : 2"
                            class="ordered-sequence-source"
                            v-on:change="replaceDataEntrySequenceNode(item.nodeIndex, $event)"
                        />
                      </div>
                    </div>
                  </div>

                </section>
              </div>
            </div>
          </template>

          <div v-else class="data-entry-file-view">
            <div class="data-entry-script-header">
              <h3>Generated input_file (.inp)</h3>
              <span>{{ generatedInputDisplayName }}</span>
            </div>
            <b-form-textarea
                v-model="inpcContent"
                class="data-entry-script-preview data-entry-file-textarea"
                rows="18"
                no-resize
                v-on:input="onFileViewInput"
            />
          </div>
        </div>
      </section>

      <RunResourceSettings
          ref="runResourceSettings"
          :run="run"
          :input-state="inputState"
          :application-module-id="resourceApplicationModuleId"
          v-on:updateResources="updateResources"
      />

      <div class="validation-copy" aria-hidden="true">
        <b-form-invalid-feedback :state="inputState.root">
          Run root is required
        </b-form-invalid-feedback>
        <b-form-invalid-feedback :state="inputState.inpcContent">
          Run input is required
        </b-form-invalid-feedback>
      </div>

      <div class="run-actions-bar">
        <button-overlay :show="processing">
          <b-button variant="outline-primary" class="run-action-button" :disabled="processing"
                    v-on:click="onSave(false)">
            Save
          </b-button>
        </button-overlay>
        <button-overlay :show="processing">
          <b-button variant="primary" class="run-action-button" :disabled="processing" v-on:click="onSave(true)">
            Submit
          </b-button>
        </button-overlay>
      </div>
    </main>

    <b-modal
        v-model="showWorkflowStepGuardModal"
        title="Workflow step incomplete"
        ok-only
        ok-title="OK"
        ok-variant="primary"
        centered
        return-focus="#new-run-top-anchor"
        v-on:hidden="onWorkflowStepGuardModalHidden"
    >
      {{ workflowStepGuardMessage }}
    </b-modal>

  </div>
</template>

<script>
import store from "@/store";
import router from "@/router";
import { eventBus } from "@/event-bus";
import ButtonOverlay from "@/components/overlay/button-overlay";
import RunResourceSettings from "@/components/blocks/RunResourceSettings";
import UserStorage from "@/components/overlay/UserStorage";
import { descriptions } from "@/fileData";
import { InputService } from "@/service/epolyscat-service";
import {
  appendEPolyScatSequenceNode,
  buildEPolyScatInputScript as buildEPolyScatInputScriptFromValues,
  EPOLYSCAT_DATA_ENTRY_SECTIONS,
  EPOLYSCAT_INPUT_SCHEMA,
  EPOLYSCAT_OUTPUT_DEFINITIONS,
  EPOLYSCAT_RECOMMENDED_VALUES as DEFAULT_DATA_ENTRY_VALUES,
  listEPolyScatSequenceNodes,
  moveEPolyScatSequenceNode,
  normalizeEPolyScatInputContents,
  parseEPolyScatInputScript as parseEPolyScatInputScriptFromContents,
  patchEPolyScatInputScript,
  removeEPolyScatSequenceNode,
  replaceEPolyScatSequenceNode,
} from "@/utils/epolyscat-input-script";
import { buildWorkflowOutputInputBinding } from "@/utils/workflow-file-linking";

export default {
  name: 'CreateRun',
  components: {ButtonOverlay, RunResourceSettings, UserStorage},
  store: store,
  data() {
    const storedShowDescriptions = this.$store.getters["settings/getPreference"]("showDescriptions");

    return {
      experimentId: null,
      viewIds: [],
      showDescriptions: storedShowDescriptions == null ? true : storedShowDescriptions,
      showWorkflowStepGuardModal: false,
      workflowStepGuardPendingScroll: false,
      workflowStepGuardMessage: "",
      workflowOutputBindingError: null,
      selectedRunType: "module",
      selectedApplicationId: "ePolyScat",
      selectedWorkflowApplicationId: "OpenMolcas",
      workflowStageSelections: {
        Data_Gen: "OpenMolcas",
        ePolyScat_Run: "ePolyScat",
        Analysis: "CnvMath",
      },
      runTypeOptions: [
        {
          id: "module",
          label: "Modules",
          title: "Modules",
          description: "Single application run",
          applications: [
            {
              id: "Gaussian16",
              label: "Gaussian16",
              description: "Input computation",
              requiredFiles: [
                {
                  name: "Gaussian_Input",
                  label: "Gaussian_Input",
                  description: "Gaussian input file",
                  generatedFileName: "gaussian.inp",
                },
              ],
            },
            {
              id: "OpenMolcas",
              label: "OpenMolcas",
              description: "Input computation",
              requiredFiles: [
                {
                  name: "Molcas_Input",
                  label: "Molcas_Input",
                  description: "Molcas input file",
                  generatedFileName: "molcas.input",
                },
              ],
            },
            {
              id: "ePolyScat",
              label: "ePolyScat",
              description: "Scattering calculation",
              requiredFiles: [
                {
                  name: "ePolyScat_Input_Data",
                  label: "input_data",
                  description: "Quantum chemistry output data",
                  generatedFileName: "input_data",
                },
                {
                  name: "ePolyscat_Input_File",
                  label: "input_file (.inp)",
                  description: "ePolyScat data records and commands",
                  generatedFileName: "input_file.inp",
                },
              ],
            },
          ],
        },
        {
          id: "utility",
          label: "Utilities",
          title: "Utilities",
          description: "Post-processing utility",
          applications: [
            {
              id: "CnvMath",
              label: "CnvMath",
              description: "ConvertToMathematica",
              requiredFiles: [
                {
                  name: "BendOrient_Output",
                  label: "BendOrient_Output",
                  description: "BendOrient output file",
                  generatedFileName: "bendorient.dat",
                },
              ],
            },
            {
              id: "CnvMatLab",
              label: "CnvMatLab",
              description: "ConvertToMatlab",
              requiredFiles: [
                {
                  name: "BendOrient_Output",
                  label: "BendOrient_Output",
                  description: "BendOrient output file",
                  generatedFileName: "bendorient.dat",
                },
              ],
            },
            {
              id: "CnvLinFull",
              label: "CnvLinFull",
              description: "Compute Ph.ion.diff.Xsec.",
              requiredFiles: [
                {
                  name: "DumpOut",
                  label: "DumpOut",
                  description: "DumpIdy output file",
                  generatedFileName: "dumpidy.dat",
                },
              ],
            },
            {
              id: "MoldenMerge",
              label: "MoldenMerge",
              description: "Merge Molden Data Files",
              requiredFiles: [
                {
                  name: "molden.dat",
                  label: "molden.dat",
                  description: "Molden data file",
                  generatedFileName: "molden.dat",
                },
              ],
            },
            {
              id: "NRFPAD",
              label: "NRFPAD",
              description: "N photon RFPAD",
              requiredFiles: [
                {
                  name: "Cross_Section_Input_File",
                  label: "Cross_Section_Input_File",
                  description: "OrientNCro output file",
                  generatedFileName: "orientncro.dat",
                },
              ],
            },
            {
              id: "Cube2igor",
              label: "Cube2igor",
              description: "G16CubeToIGOR Plots",
              requiredFiles: [
                {
                  name: "Cube_Output",
                  label: "Cube_Output",
                  description: "Gaussian cube output file",
                  generatedFileName: "gaussian.cube",
                },
              ],
            },
          ],
        },
        {
          id: "workflow",
          label: "Workflows",
          title: "Workflow Stages",
          description: "Staged calculation",
          applications: [
            {
              id: "Data_Gen",
              label: "Data Generation",
              description: "Gaussian16 or OpenMolcas",
              workflowApplicationIds: ["Gaussian16", "OpenMolcas"],
            },
            {
              id: "ePolyScat_Run",
              label: "ePolyScat Run",
              description: "Run ePolyScat",
              workflowApplicationIds: ["ePolyScat"],
            },
            {
              id: "Analysis",
              label: "Post-processing",
              description: "Analysis utilities",
              workflowApplicationIds: [
                "CnvMath",
                "CnvMatLab",
                "CnvLinFull",
                "MoldenMerge",
                "NRFPAD",
                "Cube2igor",
              ],
            },
            {
              id: "Visualization",
              label: "Visualization",
              description: "View and plot completed outputs",
              workflowApplicationIds: [],
              localOnly: true,
            },
          ],
        },
      ],
      run: {
        name: `EPolyScat Run on ${(new Date(Date.now())).toLocaleString()}`,
        id: -1,
        status: "UNSUBMITTED",
        jobStatus: "UNSUBMITTED",
        displayStatus: "UNSUBMITTED",
        isEmailNotificationOn: false,
        description: "",
        groupResourceProfileId: null,
        computeResourceId: null,
        queueName: null,
        coreCount: null,
        nodeCount: null,
        wallTimeLimit: null,
        totalPhysicalMemory: 0,
        executions: [],
        viewIds: []
      },
      selectedInputFile: "ePolyscat_Input_File",
      selectedDataEntryTab: "grid-expansion",
      dataEntryViewMode: "table",
      sequenceNodeType: "DataRecord",
      sequenceNodeLabel: "LMax",
      sequenceNodeTypeOptions: [
        { value: "DataRecord", text: "Data record" },
        { value: "Command", text: "Command" },
      ],

      inpcContent: null,
      inpcContentType: 'text',

      inputFiles: [
        {
          id: "input_data",
          label: "input_data",
          inputFileName: "ePolyScat_Input_Data",
          generatedFileName: "input_data",
        },
        {
          id: "input_file",
          label: "input_file (.inp)",
          inputFileName: "ePolyscat_Input_File",
          generatedFileName: "input_file.inp",
        },
      ],
      dataEntrySections: EPOLYSCAT_DATA_ENTRY_SECTIONS,
      outputDefinitions: EPOLYSCAT_OUTPUT_DEFINITIONS,
      dataEntryRecommendedValues: { ...DEFAULT_DATA_ENTRY_VALUES },
      dataEntryValues: { ...DEFAULT_DATA_ENTRY_VALUES },
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

      isApplyingFileContent: false,
      processing: false
    };
  },
  computed: {
    pageTitle() {
      if (this.isWorkflowChildEdit) {
        return "Workflow Step";
      }
      return this.cloneRunId ? "Clone Run" : "New Run";
    },
    workflowOutputBindingMessage() {
      const error = this.workflowOutputBindingError;
      if (!error) {
        return "";
      }
      if (typeof error === "string") {
        return error;
      }
      return error.message
          || (error.response && error.response.data && error.response.data.detail)
          || "The previous run output could not be inherited.";
    },
    isWorkflowChildEdit() {
      return Boolean(this.workflowChildRunId);
    },
    activeRunTypeOption() {
      return this.runTypeOptions.find(runType => runType.id === this.selectedRunType) || this.runTypeOptions[0];
    },
    activeRunTypeTitle() {
      return this.activeRunTypeOption ? this.activeRunTypeOption.title : "";
    },
    activeRunTypeApplications() {
      return this.activeRunTypeOption ? this.activeRunTypeOption.applications : [];
    },
    activeRunApplication() {
      return this.activeRunTypeApplications.find(application => application.id === this.selectedApplicationId)
          || this.activeRunTypeApplications[0]
          || null;
    },
    moduleApplications() {
      const moduleOption = this.runTypeOptions.find(runType => runType.id === "module");
      return moduleOption ? moduleOption.applications : [];
    },
    utilityApplications() {
      const utilityOption = this.runTypeOptions.find(runType => runType.id === "utility");
      return utilityOption ? utilityOption.applications : [];
    },
    reusableApplications() {
      return [...this.moduleApplications, ...this.utilityApplications];
    },
    activeWorkflowApplications() {
      if (this.selectedRunType !== "workflow" || !this.activeRunApplication) {
        return [];
      }

      return (this.activeRunApplication.workflowApplicationIds || [])
          .map(applicationId => this.reusableApplications.find(application => application.id === applicationId))
          .filter(application => !!application);
    },
    activeWorkflowApplicationTitle() {
      if (this.selectedApplicationId === "Data_Gen") {
        return "Input Computation";
      }

      if (this.selectedApplicationId === "Analysis") {
        return "Post-processing Utility";
      }

      return "Application";
    },
    activeExecutionApplication() {
      if (this.selectedRunType !== "workflow") {
        return this.activeRunApplication;
      }

      return this.activeWorkflowApplications.find(
          application => application.id === this.activeWorkflowApplicationId
      ) || this.activeWorkflowApplications[0] || this.activeRunApplication;
    },
    resourceApplicationModuleId() {
      const epolyscatApplicationModuleId = this.$store.getters[
          "settings/epolyscatApplicationModuleId"
      ];
      if (this.activeExecutionApplication && this.activeExecutionApplication.id === "Gaussian16") {
        return this.$store.getters["settings/gaussian16ApplicationModuleId"]
            || epolyscatApplicationModuleId;
      }
      if (this.activeExecutionApplication && this.activeExecutionApplication.id === "OpenMolcas") {
        return this.$store.getters["settings/openmolcasApplicationModuleId"]
            || epolyscatApplicationModuleId;
      }
      return epolyscatApplicationModuleId;
    },
    activeWorkflowApplicationId() {
      if (this.selectedRunType !== "workflow") {
        return "";
      }

      return this.workflowStageSelections[this.selectedApplicationId]
          || this.selectedWorkflowApplicationId
          || this.getDefaultWorkflowApplicationId();
    },
    activeRequiredFiles() {
      return this.activeExecutionApplication && this.activeExecutionApplication.requiredFiles
          ? this.activeExecutionApplication.requiredFiles
          : [];
    },
    hasRequiredFiles() {
      return this.activeRequiredFiles.length > 0;
    },
    descriptionPath() {
      if (this.selectedRunType === "module") {
        return ["MODULE", this.selectedApplicationId];
      }

      if (this.selectedRunType === "utility") {
        return ["UTILITY", this.selectedApplicationId];
      }

      return [
        "WORKFLOW",
        this.selectedApplicationId,
        this.activeExecutionApplication ? this.activeExecutionApplication.id : "",
      ].filter(Boolean);
    },
    newRunDescriptions() {
      return this.descriptionPath
          .filter((item, index, path) => item in descriptions && path.indexOf(item) === index)
          .map(item => ({
            name: item,
            text: descriptions[item],
          }));
    },
    selectedRunMode() {
      return this.selectedRunType;
    },
    selectedModuleApplication() {
      if (this.isWorkflowChildEdit && this.selectedApplicationId === "ePolyScat_Run") {
        return "ePolyScat";
      }
      return this.selectedRunType === "module" ? this.selectedApplicationId : "";
    },
    selectedWorkflowStage() {
      return this.isWorkflowChildEdit ? this.selectedApplicationId : "";
    },
    selectedWorkflowApplication() {
      if (this.selectedRunType !== "workflow") {
        return "";
      }

      if (this.isWorkflowChildEdit) {
        return this.selectedApplicationId === "Data_Gen" ? this.activeWorkflowApplicationId : "";
      }

      return this.workflowStageSelections.Data_Gen || "OpenMolcas";
    },
    selectedUtilityApplication() {
      if (this.selectedRunType === "utility") {
        return this.selectedApplicationId;
      }

      if (this.isWorkflowChildEdit) {
        return this.selectedApplicationId === "Analysis" ? this.activeWorkflowApplicationId : "";
      }

      return this.selectedRunType === "workflow"
          ? this.workflowStageSelections.Analysis || "CnvMath"
          : "";
    },
    selectedRunPath() {
      if (this.selectedRunType === "module") {
        return ["MODULE", this.selectedApplicationId];
      }

      if (this.selectedRunType === "utility") {
        return ["UTILITY", this.selectedApplicationId];
      }

      if (this.selectedApplicationId === "Data_Gen") {
        return ["MODULE", this.activeWorkflowApplicationId];
      }

      if (this.selectedApplicationId === "Analysis") {
        return ["UTILITY", this.activeWorkflowApplicationId];
      }

      return ["MODULE", "ePolyScat"];
    },
    root: {
      get() {
        return this.run.name;
      },
      set(value) {
        this.run.name = value;
      }
    },
    groupResourceProfileId: {
      get() {
        return this.run.groupResourceProfileId;
      },
      set(value) {
        this.run.groupResourceProfileId = value;
      }
    },
    computeResourceId: {
      get() {
        return this.run.computeResourceId;
      },
      set(value) {
        this.run.computeResourceId = value;
      }
    },
    queueName: {
      get() {
        return this.run.queueName;
      },
      set(value) {
        this.run.queueName = value;
      }
    },
    coreCount: {
      get() {
        return this.run.coreCount;
      },
      set(value) {
        this.run.coreCount = value;
      }
    },
    nodeCount: {
      get() {
        return this.run.nodeCount;
      },
      set(value) {
        this.run.nodeCount = value;
      }
    },
    wallTimeLimit: {
      get() {
        return this.run.wallTimeLimit;
      },
      set(value) {
        this.run.wallTimeLimit = value;
      }
    },
    totalPhysicalMemory: {
      get() {
        return this.run.totalPhysicalMemory;
      },
      set(value) {
        this.run.totalPhysicalMemory = value;
      }
    },
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
        inpcContent: !this.hasRequiredFiles ? true
            : this.inpcContentType !== "text" ? true
                : (!!this.inpcContent && this.inpcContent.length >= 1),
        groupResourceProfileId: !!this.groupResourceProfileId,
        computeResourceId: !!this.computeResourceId,
        queueName: !!this.queueName,
        coreCount: this.coreCount > 0,
        nodeCount: this.nodeCount > 0,
        wallTimeLimit: this.wallTimeLimit > 0,
        totalPhysicalMemory: `${this.totalPhysicalMemory}`.length >= 1
      }
    },
    isFormValid() {
      let _isFormValid = true;
      for (let i = 0; i < this.inputFieldsList.length; i++) {
        _isFormValid = _isFormValid && this.isValid[this.inputFieldsList[i]];
      }

      return _isFormValid;
    },
    activeInputFile() {
      return this.inputFiles.find(file => file.id === this.selectedInputFile) || this.inputFiles[0] || {
        id: "none",
        label: "No required file",
        inputFileName: "",
        generatedFileName: "",
      };
    },
    activeInputFiles() {
      if (!this.activeInputFile.inputFileName) {
        return [];
      }

      const { inputFiles } = this.$store.getters["input/getInputs"]({
        customPath: this.selectedRunPath
      });
      const inputFile = inputFiles[this.activeInputFile.inputFileName];

      return inputFile ? inputFile.files.filter(file => !file.deleted && !file.generatedByNewRun) : [];
    },
    activeInputDefinition() {
      if (!this.activeInputFile.inputFileName) {
        return null;
      }

      const { inputFiles } = this.$store.getters["input/getInputs"]({
        customPath: this.selectedRunPath
      });

      return inputFiles[this.activeInputFile.inputFileName] || null;
    },
    activeInputAllowsMultiple() {
      return this.activeInputDefinition ? this.activeInputDefinition.isMultiFileInput : false;
    },
    isEPolyScatScriptInput() {
      return this.activeInputFile.id === "input_file"
          || this.activeInputFile.inputFileName === "ePolyscat_Input_File";
    },
    dataEntryCardClasses() {
      return [
        `${this.dataEntryViewMode}-mode`,
        `input-${this.activeInputFile.id}`,
      ];
    },
    generatedInputFileName() {
      return this.activeInputFile.inputFileName;
    },
    generatedInputDisplayName() {
      return this.activeInputFile.generatedFileName;
    },
    activeDataEntrySection() {
      return this.dataEntrySections.find(section => section.id === this.selectedDataEntryTab)
          || this.dataEntrySections[0];
    },
    activeDataEntryFields() {
      return this.activeDataEntrySection && this.activeDataEntrySection.fields
          ? this.activeDataEntrySection.fields
          : [];
    },
    activeDataEntryTabLabel() {
      return this.activeDataEntrySection ? this.activeDataEntrySection.label : "";
    },
    dataEntrySequenceNodes() {
      return listEPolyScatSequenceNodes(this.inpcContent || "");
    },
    sequenceNodeLabelOptions() {
      const schemaItems = this.sequenceNodeType === "Command"
          ? EPOLYSCAT_INPUT_SCHEMA.commands
          : EPOLYSCAT_INPUT_SCHEMA.dataRecords;
      return schemaItems.map(item => ({ value: item.label, text: item.label }));
    },
    ePolyScatInputScript() {
      return this.buildEPolyScatInputScript();
    },
    newRunLink() {
      let newRunLink = "/runs/new?";
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
      return this.$route.query.cloneRunId || this.$route.query.clonedFrom;
    },
    workflowChildRunId() {
      return this.$route.query.workflowChildRunId;
    },
    workflowParentRunId() {
      return this.$route.query.workflowParentRunId;
    },
    outputsRunId() {
      return this.$route.query.withOutputsFrom;
    },
    cloneRunOriginalName() {
      const run = this.$store.getters["run/getRun"](this.cloneRunId);
      if (run) {
        return run.name
      } else {
        return null;
      }
    },
    experiments() {
      return this.$store.getters["experiment/getExperiments"]();
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
      if (this.workflowChildRunId) {
        return this.workflowChildRunId;
      } else if (this.cloneRunId) {
        return this.cloneRunId;
      } else if (this.runs && this.runs.length > 0) {
        return this.runs[this.runs.length - 1].runId;
      } else {
        return null;
      }
    },
    latestRun() {
      return this.$store.getters["run/getRun"](this.latestRunId);
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
  },
  methods: {
    selectRunType(runTypeId) {
      this.rememberWorkflowStageSelection();
      this.selectedRunType = runTypeId;

      const firstApplication = this.activeRunTypeApplications[0];
      this.selectedApplicationId = firstApplication ? firstApplication.id : "";
      this.selectedWorkflowApplicationId = this.getWorkflowStageSelection(this.selectedApplicationId);
      this.applySelectedRunConfiguration();
    },
    selectRunApplication(applicationId) {
      const application = this.activeRunTypeApplications.find(item => item.id === applicationId);
      if (application && application.localOnly) {
        return;
      }
      this.rememberWorkflowStageSelection();
      this.selectedApplicationId = applicationId;
      this.selectedWorkflowApplicationId = this.getWorkflowStageSelection(applicationId);
      this.applySelectedRunConfiguration();
    },
    selectWorkflowApplication(applicationId) {
      this.selectedWorkflowApplicationId = applicationId;
      this.workflowStageSelections = {
        ...this.workflowStageSelections,
        [this.selectedApplicationId]: applicationId,
      };
      this.applySelectedRunConfiguration();
    },
    rememberWorkflowStageSelection() {
      if (this.selectedRunType !== "workflow" || !this.selectedApplicationId) {
        return;
      }

      this.workflowStageSelections = {
        ...this.workflowStageSelections,
        [this.selectedApplicationId]: this.selectedWorkflowApplicationId || this.getDefaultWorkflowApplicationId(),
      };
    },
    getWorkflowStageSelection(applicationId) {
      return this.workflowStageSelections[applicationId] || this.getDefaultWorkflowApplicationId(applicationId);
    },
    getDefaultWorkflowApplicationId(applicationId = this.selectedApplicationId) {
      if (this.selectedRunType !== "workflow") {
        return "";
      }

      const workflowStage = this.activeRunTypeApplications.find(application => application.id === applicationId);
      const workflowApplicationIds = workflowStage ? workflowStage.workflowApplicationIds || [] : [];

      if (applicationId === "Data_Gen" && workflowApplicationIds.includes("OpenMolcas")) {
        return "OpenMolcas";
      }

      return workflowApplicationIds[0] || "";
    },
    applySelectedRunConfiguration({ replaceExisting = false } = {}) {
      this.syncRequiredInputFiles();
      this.$store.commit("input/SET_PATH", { path: this.selectedRunPath });
      this.syncGeneratedInpcContent();
      this.syncAllGeneratedInputFiles({ replaceExisting });
    },
    syncRequiredInputFiles() {
      this.inputFiles = this.activeRequiredFiles.map(file => ({
        id: file.name,
        label: file.label,
        inputFileName: file.name,
        generatedFileName: file.generatedFileName || file.name,
      }));

      if (!this.inputFiles.some(file => file.id === this.selectedInputFile)) {
        const preferredInputFile = this.inputFiles.find(file => file.inputFileName === "ePolyscat_Input_File");
        this.selectedInputFile = preferredInputFile
            ? preferredInputFile.id
            : this.inputFiles.length > 0 ? this.inputFiles[0].id : null;
      }
    },
    requiredFileAllowsMultiple(fileName) {
      const { inputFiles } = this.$store.getters["input/getInputs"]({
        customPath: this.selectedRunPath
      });
      const inputFile = inputFiles[fileName];

      return inputFile ? inputFile.isMultiFileInput : false;
    },
    uploadedRequiredFiles(fileName) {
      const { inputFiles } = this.$store.getters["input/getInputs"]({
        customPath: this.selectedRunPath
      });
      const inputFile = inputFiles[fileName];

      return inputFile ? inputFile.files.filter(file => !file.deleted && !file.generatedByNewRun) : [];
    },
    workflowStageRunPath(stageId) {
      const applicationId = this.getWorkflowStageSelection(stageId);

      if (stageId === "Data_Gen") {
        return ["MODULE", applicationId];
      }

      if (stageId === "Analysis") {
        return ["UTILITY", applicationId];
      }

      return ["MODULE", "ePolyScat"];
    },
    workflowStageRequiredFiles(stageId) {
      const applicationId = this.getWorkflowStageSelection(stageId);
      const application = this.reusableApplications.find(item => item.id === applicationId);

      return application && application.requiredFiles ? application.requiredFiles : [];
    },
    uploadedWorkflowStageRequiredFiles(stageId, fileName) {
      const { inputFiles } = this.$store.getters["input/getInputs"]({
        customPath: this.workflowStageRunPath(stageId)
      });
      const inputFile = inputFiles[fileName];

      return inputFile ? inputFile.files.filter(file => !file.deleted && !file.generatedByNewRun) : [];
    },
    firstIncompleteWorkflowStage() {
      if (this.selectedRunType !== "workflow") {
        return null;
      }

      return this.activeRunTypeApplications.find(stage => {
        const requiredFiles = this.workflowStageRequiredFiles(stage.id);

        return requiredFiles.some(file => this.uploadedWorkflowStageRequiredFiles(stage.id, file.name).length === 0);
      }) || null;
    },
    firstWorkflowSubmitStage() {
      return this.selectedRunType === "workflow" && this.activeRunTypeApplications.length > 0
          ? this.activeRunTypeApplications[0]
          : null;
    },
    validateWorkflowSubmitStep() {
      if (this.selectedRunType !== "workflow") {
        return true;
      }
      if (this.isWorkflowChildEdit) {
        return true;
      }

      const firstSubmitStage = this.firstWorkflowSubmitStage();
      if (!firstSubmitStage || this.selectedApplicationId === firstSubmitStage.id) {
        return true;
      }

      this.selectRunApplication(firstSubmitStage.id);
      this.workflowStepGuardMessage = `Please complete ${firstSubmitStage.label} before submitting this workflow.`;
      this.workflowStepGuardPendingScroll = true;
      this.showWorkflowStepGuardModal = true;
      return false;
    },
    onWorkflowStepGuardModalHidden() {
      if (!this.workflowStepGuardPendingScroll) {
        return;
      }

      this.workflowStepGuardPendingScroll = false;
      window.requestAnimationFrame(() => {
        window.setTimeout(() => {
          this.scrollToNewRunTop();
        }, 0);
      });
    },
    scrollToNewRunTop() {
      this.$nextTick(() => {
        if (document.activeElement && document.activeElement.blur) {
          document.activeElement.blur();
        }

        const scrollContainer = this.$refs.newRunContent && this.$refs.newRunContent.closest
            ? this.$refs.newRunContent.closest(".overflow-auto")
            : null;

        if (scrollContainer && scrollContainer.scrollTo) {
          scrollContainer.scrollTo({
            top: 0,
            behavior: "smooth",
          });
          return;
        }

        if (this.$refs.newRunContent && this.$refs.newRunContent.scrollIntoView) {
          this.$refs.newRunContent.scrollIntoView({
            behavior: "smooth",
            block: "start",
          });
        }
      });
    },
    selectRequiredFile(fileName) {
      this.selectInputFile(fileName);
      this.scrollToInputFiles();
    },
    scrollToInputFiles() {
      this.$nextTick(() => {
        if (this.$refs.inputFilesPanel) {
          this.$refs.inputFilesPanel.scrollIntoView({
            behavior: "smooth",
            block: "start",
          });
        }
      });
    },
    async onInputStorageFilesSelected(files) {
      const selectedFiles = (files || []).map(file => ({
        ...file,
        isFromComputer: false,
      }));

      await this.addFilesToActiveInput(selectedFiles);
    },
    async onInputComputerFilesSelected(event) {
      const selectedFiles = Array.from(event.target.files || []).map(file => {
        file.isFromComputer = true;
        return file;
      });

      await this.addFilesToActiveInput(selectedFiles);
      event.target.value = "";
    },
    async addFilesToActiveInput(files) {
      const inputFileName = this.activeInputFile.inputFileName;
      if (!inputFileName || files.length === 0) {
        return;
      }

      const allowsMultiple = this.activeInputAllowsMultiple;
      const filesToAdd = allowsMultiple ? files : files.slice(0, 1);
      const existingFiles = this.uploadedRequiredFiles(inputFileName);

      const preparedFiles = await Promise.all(filesToAdd.map(async (file, index) => {
        const preparedFile = file.isFromComputer ? file : { ...file };

        Object.assign(preparedFile, {
          replaceCurrent: !allowsMultiple && (existingFiles.length > 0 || index > 0),
          deleted: false,
          unchanged: false,
        });

        if (inputFileName === "ePolyscat_Input_File") {
          const contents = await this.readInputFileContents(file);
          if (contents != null) {
            preparedFile.contents = normalizeEPolyScatInputContents(contents);
          }
        }

        return preparedFile;
      }));

      preparedFiles.forEach(file => {
        this.$store.commit("input/ADD_TO_INPUT_FILE", {
          inputFileName,
          file
        });
      });

      if (inputFileName === "ePolyscat_Input_File") {
        await this.loadInpcContentFromFile(preparedFiles[0]);
      }
    },
    addWorkflowOutputToInput(inputFileName, outputFile) {
      if (!inputFileName || !outputFile) {
        return;
      }

      const existingFiles = this.uploadedRequiredFiles(inputFileName);
      this.$store.commit("input/ADD_TO_INPUT_FILE", {
        inputFileName,
        file: {
          ...outputFile,
          dataProductURI: outputFile.dataProductURI || outputFile["data-product-uri"],
          deleted: false,
          unchanged: false,
          replaceCurrent: existingFiles.length > 0,
        },
      });
    },
    workflowApplicationIdFromRun(run) {
      return run
          ? run.workflowApplication || run.moduleApplication || run.utilityApplication || ""
          : "";
    },
    applyWorkflowDataEntryValues(dataEntryValues) {
      if (!dataEntryValues) {
        return;
      }

      Object.assign(this.dataEntryValues, dataEntryValues);
      this.patchInpcContentFromTableValues(Object.keys(dataEntryValues));
    },
    async loadWorkflowPreviousOutputsIntoActiveInputs(outputsRunId, sourceRun = null) {
      this.workflowOutputBindingError = null;
      try {
        if (this.selectedApplicationId === "ePolyScat_Run") {
          const backendBinding = await InputService.fetchWorkflowOutputBinding({
            runId: outputsRunId,
            targetStageId: this.selectedApplicationId,
            targetApplicationId: this.activeWorkflowApplicationId,
            requiredFileName: "ePolyScat_Input_Data",
          });
          if (backendBinding.status === "blocked") {
            this.workflowOutputBindingError = backendBinding.message;
            return;
          }
          if (backendBinding.status === "ready") {
            this.addWorkflowOutputToInput(
                backendBinding.inputFileName,
                backendBinding.outputFile,
            );
            this.applyWorkflowDataEntryValues(backendBinding.dataEntryValues);
          }
          return;
        }

        for (const requiredFile of this.activeRequiredFiles) {
          const requiredFileName = requiredFile.name;
          const backendBinding = await InputService.fetchWorkflowOutputBinding({
            runId: outputsRunId,
            targetStageId: this.selectedApplicationId,
            targetApplicationId: this.activeWorkflowApplicationId,
            requiredFileName,
          });
          if (backendBinding.status === "blocked") {
            this.workflowOutputBindingError = backendBinding.message;
            return;
          }
          if (backendBinding.status === "ready") {
            this.addWorkflowOutputToInput(
                backendBinding.inputFileName,
                backendBinding.outputFile,
            );
            this.applyWorkflowDataEntryValues(backendBinding.dataEntryValues);
          }
        }
        return;
      } catch (error) {
        this.workflowOutputBindingError = error;
      }

      const outputFiles = await InputService.fetchOutputs(outputsRunId);
      const byName = outputFiles.reduce((filesByName, file) => ({
        ...filesByName,
        [file.name]: file,
      }), {});
      const sourceApplicationId = this.workflowApplicationIdFromRun(sourceRun);

      if (this.selectedApplicationId === "ePolyScat_Run") {
        const binding = buildWorkflowOutputInputBinding({
          outputFiles,
          targetStageId: this.selectedApplicationId,
          sourceApplicationId,
        });

        if (binding) {
          this.addWorkflowOutputToInput(binding.inputFileName, binding.outputFile);
          this.applyWorkflowDataEntryValues(binding.dataEntryValues);
        }
        return;
      }

      this.activeRequiredFiles.forEach(requiredFile => {
        const binding = buildWorkflowOutputInputBinding({
          outputFiles,
          targetStageId: this.selectedApplicationId,
          targetApplicationId: this.activeWorkflowApplicationId,
          sourceApplicationId,
          requiredFileName: requiredFile.name,
        });

        if (binding) {
          this.addWorkflowOutputToInput(binding.inputFileName, binding.outputFile);
          this.applyWorkflowDataEntryValues(binding.dataEntryValues);
        } else if (byName[requiredFile.name]) {
          this.addWorkflowOutputToInput(requiredFile.name, byName[requiredFile.name]);
        }
      });
    },
    removeInputFile(file) {
      if (!file || !this.activeInputFile.inputFileName) {
        return;
      }

      const inputFileName = this.activeInputFile.inputFileName;
      this.$store.commit("input/REMOVE_FILE", {
        filename: file.name,
        inputFileName: this.activeInputFile.inputFileName,
      });

      if (inputFileName === "ePolyscat_Input_File") {
        this.resetDataEntryValuesToDefaults();
        this.inpcContent = this.buildEPolyScatInputScript();
        this.syncInpcContentToActiveInputFile();
      }
    },
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
    async selectInputFile(fileId) {
      this.selectedInputFile = fileId;
      if (this.isEPolyScatScriptInput) {
        await this.loadInpcContentFromActiveInputFile();
      } else {
        this.syncGeneratedInpcContent();
      }
    },
    selectDataEntryViewMode(mode) {
      if (mode === "table") {
        this.applyInpcContentToDataEntryValues(this.inpcContent, { releaseLock: false });
        this.dataEntryViewMode = mode;
        this.releaseFileContentApplyLock();
        return;
      }

      this.dataEntryViewMode = mode;
    },
    selectDataEntryTab(tabId) {
      this.selectedDataEntryTab = tabId;
    },
    onSequenceNodeTypeInput(type) {
      this.sequenceNodeType = type;
      const firstOption = this.sequenceNodeLabelOptions[0];
      this.sequenceNodeLabel = firstOption ? firstOption.value : "";
    },
    sequenceNodeTypeText(type) {
      if (type === "DataRecord") return "Data record";
      if (type === "Command") return "Command";
      return "Unrecognized input";
    },
    applyDataEntrySequenceContents(contents) {
      this.inpcContentType = "text";
      this.inpcContent = contents;
      this.applyInpcContentToDataEntryValues(contents);
      this.syncInpcContentToActiveInputFile();
    },
    appendDataEntrySequenceNode() {
      const contents = appendEPolyScatSequenceNode(this.inpcContent || "", {
        type: this.sequenceNodeType,
        label: this.sequenceNodeLabel,
      });
      this.applyDataEntrySequenceContents(contents);
    },
    replaceDataEntrySequenceNode(nodeIndex, raw) {
      const contents = replaceEPolyScatSequenceNode(this.inpcContent || "", nodeIndex, raw);
      this.applyDataEntrySequenceContents(contents);
    },
    removeDataEntrySequenceNode(nodeIndex) {
      const contents = removeEPolyScatSequenceNode(this.inpcContent || "", nodeIndex);
      this.applyDataEntrySequenceContents(contents);
    },
    moveDataEntrySequenceNode(nodeIndex, direction) {
      const contents = moveEPolyScatSequenceNode(this.inpcContent || "", nodeIndex, direction);
      this.applyDataEntrySequenceContents(contents);
    },
    fieldSelectOptions(field) {
      return [
        {
          value: "",
          text: this.dataEntryPlaceholder(field.key) || "Select",
          disabled: true,
        },
        ...(field.options || []),
      ];
    },
    dataEntryPlaceholder(fieldKey) {
      const value = this.dataEntryRecommendedValues[fieldKey];
      return value === undefined || value === null || value === "" ? "" : `Recommended: ${value}`;
    },
    blankDataEntryValues() {
      return Object.keys(this.dataEntryRecommendedValues).reduce((values, key) => {
        values[key] = "";
        return values;
      }, {});
    },
    resetDataEntryValuesToDefaults() {
      this.dataEntryValues = { ...this.dataEntryRecommendedValues };
    },
    onFileViewInput(value) {
      this.inpcContentType = "text";
      this.inpcContent = normalizeEPolyScatInputContents(value);
      this.applyInpcContentToDataEntryValues(this.inpcContent);
      this.syncInpcContentToActiveInputFile();
    },
    onTableFieldInput(fieldKey, value) {
      if (fieldKey) {
        this.dataEntryValues[fieldKey] = value;
      }

      if (this.isApplyingFileContent) {
        return;
      }

      this.patchInpcContentFromTableValues(fieldKey ? [fieldKey] : []);
    },
    patchInpcContentFromTableValues(changedKeys = []) {
      if (!this.isEPolyScatScriptInput) {
        this.syncGeneratedInpcContent();
        return;
      }

      this.inpcContentType = "text";
      this.inpcContent = patchEPolyScatInputScript(
          this.inpcContent || this.buildInpcContent(),
          this.dataEntryValues,
          this.outputDefinitions,
          { changedKeys }
      );
      this.syncInpcContentToActiveInputFile();
    },
    updateResources(resources) {
      Object.assign(this.run, resources);
    },
    async readInputFileContents(file) {
      if (!file) {
        return null;
      }

      const dataProductURI = file.dataProductURI || file["data-product-uri"];
      const readableFile = file.isFromComputer ? file : {
        ...file,
        dataProductURI,
        downloadURL: file.downloadURL || (dataProductURI
            ? `sdk/download/?data-product-uri=${encodeURIComponent(dataProductURI)}`
            : file.downloadURL),
      };

      try {
        return await InputService.fetchFileContents(readableFile);
      } catch (error) {
        eventBus.$emit("error", { name: `Error while trying to read "${file.name}"`, error });
        return null;
      }
    },
    async loadInpcContentFromActiveInputFile() {
      const file = this.getActiveInpcStoreFile();
      if (file) {
        await this.loadInpcContentFromFile(file);
      } else {
        this.syncGeneratedInpcContent();
      }
    },
    async loadInpcContentFromFile(file) {
      const contents = file && "contents" in file
          ? file.contents
          : await this.readInputFileContents(file);

      if (contents == null) {
        return;
      }

      const normalizedContents = normalizeEPolyScatInputContents(contents);
      this.inpcContentType = "text";
      this.inpcContent = normalizedContents;
      this.applyInpcContentToDataEntryValues(normalizedContents);
      this.syncInpcContentToActiveInputFile();
    },
    getActiveInpcStoreFile() {
      if (!this.isEPolyScatScriptInput || !this.activeInputFile.inputFileName) {
        return null;
      }

      const { inputFiles } = this.$store.getters["input/getInputs"]({
        customPath: this.selectedRunPath,
        removeDeleted: false,
      });
      const inputFile = inputFiles[this.activeInputFile.inputFileName];

      if (!inputFile) {
        return null;
      }

      return inputFile.files.find(file => !file.deleted && !file.generatedByNewRun)
          || inputFile.files.find(file => !file.deleted && file.generatedByNewRun)
          || null;
    },
    applyInpcContentToDataEntryValues(contents, { releaseLock = true } = {}) {
      if (contents == null) {
        return;
      }

      const parsedValues = this.parseEPolyScatInputScript(contents);
      this.isApplyingFileContent = true;
      Object.assign(this.dataEntryValues, this.blankDataEntryValues(), parsedValues);
      if (releaseLock) {
        this.releaseFileContentApplyLock();
      }
    },
    releaseFileContentApplyLock() {
      this.$nextTick(() => {
        window.setTimeout(() => {
          this.isApplyingFileContent = false;
        }, 0);
      });
    },
    parseEPolyScatInputScript(contents) {
      return parseEPolyScatInputScriptFromContents(
          contents,
          this.outputDefinitions
      );
    },
    buildInpcContent() {
      return this.buildEPolyScatInputScript();
    },
    buildEPolyScatInputScript() {
      return buildEPolyScatInputScriptFromValues(this.dataEntryValues, this.outputDefinitions);
    },
    syncGeneratedInpcContent() {
      this.inpcContentType = "text";
      this.inpcContent = this.buildInpcContent();
      this.syncInpcContentToActiveInputFile();
    },
    syncAllGeneratedInputFiles({ replaceExisting = true } = {}) {
      const selectedInputFile = this.selectedInputFile;

      this.inputFiles.forEach(inputFile => {
        this.selectedInputFile = inputFile.id;
        this.inpcContent = this.buildInpcContent();
        this.syncGeneratedInputFile({ replaceExisting });
      });

      this.selectedInputFile = selectedInputFile;
      this.inpcContent = this.buildInpcContent();
    },
    syncGeneratedInputFile({ replaceExisting = true } = {}) {
      this.syncInpcContentToActiveInputFile({ replaceExisting });
    },
    syncInpcContentToActiveInputFile({ replaceExisting = true } = {}) {
      const { inputFiles } = this.$store.getters["input/getInputs"]({
        customPath: this.selectedRunPath
      });
      const inputFile = inputFiles[this.generatedInputFileName];

      if (!inputFile || this.inpcContent == null || !this.isEPolyScatScriptInput) return;

      const existingFile = inputFile.files.find(file =>
          !file.deleted &&
          !file.generatedByNewRun &&
          (file.isFromComputer || file.dataProductURI || file.name !== this.generatedInputDisplayName)
      ) || inputFile.files.find(file => !file.deleted && file.generatedByNewRun);

      const hasExistingFile = inputFile.files.some(file => !file.deleted);
      const file = {
        ...(existingFile || {}),
        name: existingFile ? existingFile.name : this.generatedInputDisplayName,
        contents: this.inpcContent,
        isFromComputer: existingFile ? !!existingFile.isFromComputer : false,
        generatedByNewRun: existingFile ? !!existingFile.generatedByNewRun : true,
        deleted: false,
        unchanged: false,
        replaceCurrent: replaceExisting && hasExistingFile,
      };

      this.$store.commit("input/ADD_TO_INPUT_FILE", {
        inputFileName: this.generatedInputFileName,
        file
      });
    },
    async onSave(submit = false) {
      if (submit) {
        if (!this.validateWorkflowSubmitStep()) return;
      }

      if (this.isEPolyScatScriptInput) {
        if (this.dataEntryViewMode === "file") {
          this.applyInpcContentToDataEntryValues(this.inpcContent);
        }
        this.syncInpcContentToActiveInputFile();
      } else {
        this.syncGeneratedInpcContent();
        this.syncAllGeneratedInputFiles();
      }
      this.makeFormVisited();

      if (!this.isFormValid) return;

      this.processing = true;

      try {
        const runPayload = {
          ...this.run,
          runMode: this.selectedRunMode,
          moduleApplication: this.selectedModuleApplication,
          workflowStage: this.selectedWorkflowStage,
          workflowApplication: this.selectedWorkflowApplication,
          utilityApplication: this.selectedUtilityApplication,
          workflowMetadata: {
            isWorkflowPlan: this.selectedRunType === "workflow",
            selectedApplication: this.activeExecutionApplication ? this.activeExecutionApplication.id : "",
            requiredFiles: this.activeRequiredFiles.map(file => file.name),
            dataGenerationApplication: this.workflowStageSelections.Data_Gen || "OpenMolcas",
            analysisApplications: [this.workflowStageSelections.Analysis || "CnvMath"],
            plannedStageIds: ["Data_Gen", "ePolyScat_Run", "Analysis", "Visualization"],
          },
          viewIds: this.viewIds
        };

        if (this.isWorkflowChildEdit) {
          runPayload.workflowMetadata = {
            ...(this.run.workflowMetadata || {}),
            parent_run_id: this.run.parentRunId || this.workflowParentRunId,
            stage: this.selectedWorkflowStage,
            application: this.activeExecutionApplication ? this.activeExecutionApplication.id : "",
            requiredFiles: this.activeRequiredFiles.map(file => file.name),
          };
          this.run = await this.$store.dispatch("run/updateRun", runPayload);
        } else {
          this.run = await this.$store.dispatch("run/createRun", runPayload);
        }

        if (submit) {
          this.run = await this.$store.dispatch("run/submitRun", {
            runId: this.run.id
          });
        } else {
          this.$store.dispatch("view/fetchViews");
        }

        this.$store.dispatch("experiment/fetchExperiments");
        router.push(`/runs/${this.run.parentRunId || this.workflowParentRunId || this.run.id}`);
      } catch (error) {
        eventBus.$emit("error", { name: `Error while trying to ${submit ? "submit" : "create"} the run`, error });
      } finally {
        this.processing = false;
      }
    },
    async initializeInputBinding() {
      await this.$store.dispatch("input/fetchPathLabels");
      this.applySelectedRunConfiguration({ replaceExisting: false });
      await this.fetchURLData();
      this.syncAllGeneratedInputFiles({ replaceExisting: false });
    },
    applyWorkflowChildRun(sourceRun) {
      const stageId = sourceRun.workflowStage || "ePolyScat_Run";
      const applicationId = sourceRun.workflowApplication
          || sourceRun.moduleApplication
          || sourceRun.utilityApplication
          || this.getDefaultWorkflowApplicationId(stageId);

      this.selectedRunType = "workflow";
      this.selectedApplicationId = stageId;
      this.selectedWorkflowApplicationId = applicationId;
      this.workflowStageSelections = {
        ...this.workflowStageSelections,
        [stageId]: applicationId,
      };
      this.run = {
        ...this.run,
        ...sourceRun,
      };
      this.viewIds = sourceRun.isTutorial ? [] : sourceRun.viewIds || [];
      this.applySelectedRunConfiguration({ replaceExisting: false });
    },
    async fetchURLData() {
      if (this.viewId) {
        this.viewIds = [this.viewId];
      }

      if (this.workflowChildRunId) {
        const workflowChildRunId = parseInt(this.workflowChildRunId);
        const sourceRun = await this.$store.dispatch("run/fetchRun", { runId: workflowChildRunId });

        this.applyWorkflowChildRun(sourceRun);
        await this.$store.dispatch("run/loadInputs", { runId: workflowChildRunId });

        if (this.outputsRunId) {
          const outputsRunId = parseInt(this.outputsRunId);
          const outputSourceRun = await this.$store.dispatch("run/fetchRun", { runId: outputsRunId });
          await this.loadWorkflowPreviousOutputsIntoActiveInputs(outputsRunId, outputSourceRun);
        }
        this.$store.commit("input/SET_PATH", { path: this.selectedRunPath });
      } else if (this.outputsRunId) {
        const outputsRunId = parseInt(this.outputsRunId);
        const sourceRun = await this.$store.dispatch("run/fetchRun", { runId: outputsRunId });

        await this.$store.dispatch("input/loadOutputsIntoInputs", { runId: outputsRunId });
        this.run.name = `Continuation of ${sourceRun.name}`;
        this.viewIds = sourceRun.isTutorial ? [] : sourceRun.viewIds || [];
      } else if (this.cloneRunId) {
        const clonedRunId = parseInt(this.cloneRunId);
        const sourceRun = await this.$store.dispatch("run/fetchRun", { runId: clonedRunId });

        await this.$store.dispatch("run/loadInputs", { runId: clonedRunId });
        this.viewIds = sourceRun.isTutorial ? [] : sourceRun.viewIds || [];
        this.run = {
          ...this.run,
          ...this.$store.getters["run/getStrippedRun"](sourceRun),
          name: `Clone of ${sourceRun.name}`,
          id: -1,
          status: "UNSUBMITTED",
          jobStatus: "UNSUBMITTED",
          displayStatus: "UNSUBMITTED",
          executions: []
        };
      }
    },
    async refreshData() {
      this.$store.dispatch("experiment/fetchExperiments");

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
        const run = await this.$store.dispatch("run/fetchRun", {runId: this.latestRunId});
        this.experimentId = run.experimentId;
      }
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
    async latestRun() {
      // Copy compute resource details if they exist
      if (this.latestRun && this.latestRun.groupResourceProfileId && this.latestRun.computeResourceId && this.latestRun.queueName) {
        this.groupResourceProfileId = this.latestRun.groupResourceProfileId;
        this.computeResourceId = this.latestRun.computeResourceId;
        this.queueName = this.latestRun.queueName;
        this.coreCount = this.latestRun.coreCount;
        this.nodeCount = this.latestRun.nodeCount;
        this.wallTimeLimit = this.latestRun.wallTimeLimit;
        this.totalPhysicalMemory = this.latestRun.totalPhysicalMemory;
        this.$nextTick(() => {
          if (this.$refs.runResourceSettings) {
            this.$refs.runResourceSettings.refreshResourceSelection();
          }
        });
      }
    },
    showDescriptions(value) {
      this.$store.commit("settings/setPreference", {
        key: "showDescriptions",
        value
      });
    }
  },
  async mounted() {
    this.experimentId = this.experimentIdFromQueryString;
    this.syncGeneratedInpcContent();

    await this.$store.dispatch("settings/fetchSettings");
    await this.initializeInputBinding();
    await this.refreshData();
  }
};
</script>

<style scoped>
.new-run-page {
  background: #ffffff;
  color: #000000;
  min-height: 100%;
}

.new-run-content {
  max-width: 1120px;
  padding: 18px 32px 90px 44px;
}

.new-run-top-anchor {
  height: 1px;
  outline: none;
  overflow: hidden;
  position: absolute;
  width: 1px;
}

.new-run-breadcrumb {
  font-size: 16px;
  margin-bottom: 16px;
}

.new-run-heading {
  align-items: center;
  display: flex;
  justify-content: space-between;
  margin: 0 0 36px;
  max-width: 1017px;
}

.new-run-title {
  font-size: 28px;
  font-weight: 700;
  line-height: 1.2;
  margin: 0;
}

.run-type-selector-section {
  margin-bottom: 34px;
  max-width: 1017px;
}

.run-type-switcher-row {
  display: grid;
  gap: 10px;
  grid-template-columns: repeat(3, minmax(0, 1fr));
  margin-bottom: 14px;
}

.run-type-switch-option {
  background: #ffffff;
  border: 1px solid #d4d9df;
  border-radius: 6px;
  color: #1f2933;
  min-height: 56px;
  padding: 10px 14px;
  text-align: left;
}

.run-type-switch-option.active {
  background: #226597;
  border-color: #226597;
  color: #ffffff;
}

.run-type-switch-label,
.run-type-switch-description {
  display: block;
  overflow-wrap: anywhere;
}

.run-type-switch-label {
  font-size: 15px;
  font-weight: 700;
  line-height: 1.2;
}

.run-type-switch-description {
  font-size: 12px;
  font-weight: 500;
  margin-top: 4px;
  opacity: 0.8;
}

.run-type-selection-grid {
  align-items: stretch;
  display: grid;
  gap: 18px;
  grid-template-columns: minmax(168px, 0.9fr) minmax(280px, 1.35fr) minmax(250px, 1.15fr);
}

.workflow-selection-grid {
  align-items: stretch;
  display: grid;
  gap: 18px;
  grid-template-columns: minmax(0, 1.5fr) minmax(250px, 0.85fr);
}

.run-selection-column {
  border: 1px solid #d8dde5;
  border-radius: 6px;
  min-width: 0;
  padding: 16px;
}

.run-selection-column h2 {
  color: #1b2b38;
  font-size: 15px;
  font-weight: 700;
  line-height: 1.2;
  margin: 0 0 12px;
}

.run-selection-option,
.workflow-stage-application-option {
  background: #f7f8fa;
  border: 1px solid #d4d9df;
  border-radius: 6px;
  color: #1f2933;
  display: block;
  font-size: 15px;
  line-height: 1.25;
  margin: 0 0 10px;
  min-height: 44px;
  padding: 10px 12px;
  text-align: left;
  width: 100%;
}

.run-selection-option.active,
.workflow-stage-application-option.active {
  background: #226597;
  border-color: #226597;
  color: #ffffff;
}

.run-selection-label {
  display: block;
  font-weight: 700;
  overflow-wrap: anywhere;
}

.run-selection-description {
  display: block;
  font-size: 12px;
  font-weight: 500;
  margin-top: 4px;
  opacity: 0.8;
  overflow-wrap: anywhere;
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

.workflow-stage-application-panel {
  background: #f7f8fa;
  border: 1px solid #d4d9df;
  border-radius: 6px;
  margin-top: 14px;
  padding: 12px;
}

.workflow-stage-application-panel h3 {
  color: #36454f;
  font-size: 13px;
  font-weight: 700;
  line-height: 1.2;
  margin: 0 0 10px;
}

.workflow-stage-application-options {
  display: grid;
  gap: 8px;
  grid-template-columns: repeat(auto-fit, minmax(128px, 1fr));
}

.workflow-stage-application-option {
  background: #ffffff;
  display: flex;
  flex-direction: column;
  justify-content: center;
  margin: 0;
  min-height: 58px;
  padding: 9px 10px;
}

.workflow-stage-application-label {
  font-size: 14px;
  font-weight: 700;
  overflow-wrap: anywhere;
}

.workflow-stage-application-description {
  font-size: 11px;
  font-weight: 500;
  margin-top: 3px;
  opacity: 0.78;
  overflow-wrap: anywhere;
}

.required-file-list {
  display: grid;
  gap: 10px;
}

.required-file-row {
  background: #f7f8fa;
  border: 1px solid #d4d9df;
  border-radius: 6px;
  cursor: pointer;
  min-width: 0;
  padding: 10px 12px;
}

.required-file-row:hover,
.required-file-row:focus {
  border-color: #226597;
  outline: 0;
}

.required-file-name {
  color: #1f2933;
  display: block;
  font-size: 14px;
  font-weight: 700;
  line-height: 1.25;
  overflow-wrap: anywhere;
  word-break: break-word;
}

.required-file-description,
.required-file-empty {
  color: #5f6f7a;
  display: block;
  font-size: 12px;
  font-weight: 500;
  line-height: 1.35;
  margin-top: 4px;
}

.required-file-status {
  border-radius: 999px;
  display: inline-flex;
  font-size: 12px;
  font-weight: 700;
  line-height: 1.2;
  margin-top: 10px;
  padding: 5px 8px;
}

.required-file-status-ready {
  background: #e8f5ee;
  color: #17633b;
}

.required-file-status-missing {
  background: #fff3e5;
  color: #9a5a00;
}

.required-file-selected-list {
  color: #53646f;
  display: block;
  font-size: 12px;
  line-height: 1.35;
  margin-top: 8px;
  overflow-wrap: anywhere;
}

.new-run-descriptions {
  border-top: 1px solid #d8dde5;
  margin: 0 0 34px;
  max-width: 1017px;
  padding-top: 14px;
}

.workflow-output-binding-warning {
  margin: -12px 0 28px;
  max-width: 1017px;
}

.new-run-description-header {
  align-items: center;
  display: flex;
  gap: 8px;
  margin-bottom: 8px;
}

.new-run-description-header h2 {
  color: #1b2b38;
  font-size: 15px;
  font-weight: 700;
  line-height: 1.2;
  margin: 0;
}

.new-run-description-toggle {
  font-size: 14px;
  line-height: 1.2;
  padding: 0;
}

.new-run-description-list {
  border-top: 1px solid #d8dde5;
  padding-top: 10px;
}

.new-run-description-row {
  color: #27343d;
  font-size: 14px;
  line-height: 1.45;
  margin-bottom: 6px;
  overflow-wrap: anywhere;
}

.input-files-panel {
  margin-bottom: 38px;
}

.input-file-actions {
  align-items: center;
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
  margin-top: 14px;
}

.input-file-storage-button,
.input-file-computer-button {
  border-radius: 6px;
  font-size: 13px;
  font-weight: 700;
  line-height: 1.2;
  margin: 0;
  min-height: 34px;
  padding: 9px 12px;
}

.input-file-computer-button {
  align-items: center;
  background: #ffffff;
  border: 1px solid #226597;
  color: #226597;
  cursor: pointer;
  display: inline-flex;
}

.input-file-or {
  color: #6a7882;
  font-size: 12px;
  font-weight: 700;
}

.input-file-computer-input {
  display: none;
}

.input-file-selected-list {
  display: flex;
  flex-wrap: wrap;
  gap: 6px;
  margin-top: 12px;
}

.input-file-selected {
  align-items: center;
  background: #e9f2f8;
  border: 1px solid #c7dce9;
  border-radius: 999px;
  color: #1f4e68;
  display: inline-flex;
  font-size: 12px;
  font-weight: 700;
  line-height: 1.2;
  max-width: 100%;
  min-height: 26px;
  overflow-wrap: anywhere;
  padding: 5px 26px 5px 8px;
  position: relative;
}

.input-file-selected-name {
  overflow-wrap: anywhere;
}

.input-file-remove {
  align-items: center;
  background: #ffffff;
  border: 1px solid #c7dce9;
  border-radius: 50%;
  color: #1f4e68;
  display: inline-flex;
  font-size: 14px;
  height: 18px;
  justify-content: center;
  line-height: 1;
  padding: 0;
  position: absolute;
  right: -6px;
  top: -7px;
  width: 18px;
}

.input-file-remove:hover,
.input-file-remove:focus {
  border-color: #226597;
  color: #226597;
  outline: 0;
}

.input-files-panel h2,
.data-entry-section h2 {
  font-size: 17px;
  font-weight: 700;
  line-height: 1.3;
  margin: 0 0 12px;
}

.input-file-tabs {
  border: 1px solid #dddddd;
  border-radius: 6px;
  display: inline-flex;
  gap: 8px;
  padding: 16px;
}

.input-file-tab,
.data-entry-tab,
.file-chip {
  background: #eeeeee;
  border: 0;
  border-radius: 6px;
  color: #000000;
  font-size: 16px;
  font-weight: 700;
  line-height: 1.2;
  min-height: 36px;
  padding: 8px 18px;
}

.input-file-tab.active,
.data-entry-tab.active,
.file-chip.active {
  background: #226597;
  color: #ffffff;
}

.section-title-row {
  align-items: center;
  display: flex;
  justify-content: space-between;
  margin-bottom: 10px;
}

.view-mode-switcher {
  display: flex;
  gap: 10px;
}

.view-mode-button {
  align-items: center;
  background: #eeeeee;
  border: 0;
  border-radius: 8px;
  color: #000000;
  display: inline-flex;
  font-size: 23px;
  height: 38px;
  justify-content: center;
  width: 38px;
}

.view-mode-button.active {
  background: #226597;
  color: #ffffff;
}

.data-entry-card {
  border: 1px solid #dddddd;
  border-radius: 6px;
  overflow: hidden;
}

.data-entry-card.table-mode {
  border-color: #d8d8d8;
}

.data-entry-tabs {
  align-items: center;
  border-bottom: 1px solid #dddddd;
  display: flex;
  flex-wrap: wrap;
  gap: 16px;
  min-height: 68px;
  padding: 16px;
}

.data-entry-body {
  min-height: 214px;
  padding: 24px 16px;
}

.data-entry-card.table-mode .data-entry-tabs {
  gap: 8px;
  min-height: 54px;
  padding: 10px 16px;
}

.data-entry-card.table-mode .data-entry-tab {
  min-height: 36px;
  padding: 8px 12px;
}

.data-entry-body-table {
  padding: 18px 16px;
}

.data-entry-editor-grid {
  display: block;
}

.data-entry-editor-panel {
  min-width: 0;
}

.data-entry-file-view {
  padding: 18px 16px;
}

.data-entry-section-heading,
.data-entry-script-header {
  align-items: baseline;
  display: flex;
  gap: 10px;
  justify-content: space-between;
  margin-bottom: 14px;
}

.data-entry-section-heading h3,
.data-entry-script-header h3 {
  color: #17212b;
  font-size: 15px;
  font-weight: 700;
  line-height: 1.25;
  margin: 0;
}

.data-entry-section-heading span,
.data-entry-script-header span {
  color: #61707c;
  font-size: 12px;
  font-weight: 700;
  line-height: 1.2;
  overflow-wrap: anywhere;
  text-align: right;
}

.manual-field-rows {
  display: grid;
  gap: 12px;
}

.manual-field-row {
  align-items: center;
  display: grid;
  gap: 12px;
  grid-template-columns: minmax(124px, 0.72fr) minmax(78px, 0.42fr) minmax(150px, 1fr);
  min-width: 0;
}

.data-entry-card.table-mode .manual-field-row {
  border-bottom: 1px solid #eeeeee;
  min-height: 54px;
  padding: 8px 0;
}

.data-entry-card.table-mode .manual-field-row:last-child {
  border-bottom: 0;
}

.manual-field-row label,
.output-record-type {
  color: #1f2933;
  font-size: 14px;
  font-weight: 700;
  line-height: 1.25;
  margin: 0;
  overflow-wrap: anywhere;
}

.manual-field-record {
  background: #eef3f7;
  border: 1px solid #d4dde5;
  border-radius: 999px;
  color: #314554;
  font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", monospace;
  font-size: 12px;
  font-weight: 700;
  line-height: 1.2;
  max-width: 100%;
  overflow-wrap: anywhere;
  padding: 5px 8px;
  text-align: center;
}

.manual-field-input {
  border-color: #b7c0c9;
  border-radius: 6px;
  color: #1f2933;
  font-size: 14px;
  font-weight: 700;
  height: 38px;
  min-width: 0;
}

.output-record-grid {
  display: grid;
  gap: 12px;
}

.output-record-row {
  align-items: center;
  display: grid;
  gap: 10px;
  grid-template-columns: minmax(112px, 0.62fr) minmax(130px, 0.82fr) minmax(130px, 1fr) 48px;
  min-width: 0;
}

.output-record-description,
.output-record-extension {
  color: #61707c;
  font-size: 12px;
  font-weight: 600;
  line-height: 1.3;
  overflow-wrap: anywhere;
}

.output-record-extension {
  text-align: right;
}

.ordered-sequence-editor {
  min-width: 0;
}

.ordered-sequence-toolbar {
  align-items: center;
  display: grid;
  gap: 10px;
  grid-template-columns: minmax(132px, 0.45fr) minmax(180px, 1fr) auto;
  margin-bottom: 16px;
}

.ordered-sequence-select {
  border-color: #b7c0c9;
  border-radius: 6px;
  font-size: 14px;
  font-weight: 700;
  height: 38px;
  min-width: 0;
}

.ordered-sequence-add {
  align-items: center;
  display: inline-flex;
  font-size: 13px;
  font-weight: 700;
  gap: 5px;
  height: 38px;
  justify-content: center;
  min-width: 82px;
}

.ordered-sequence-list {
  border-top: 1px solid #d8dde5;
}

.ordered-sequence-row {
  border-bottom: 1px solid #e3e7eb;
  padding: 12px 0;
}

.ordered-sequence-row-header {
  align-items: center;
  display: grid;
  gap: 10px;
  grid-template-columns: 28px minmax(0, 1fr) auto;
  margin-bottom: 8px;
}

.ordered-sequence-index {
  align-items: center;
  background: #eef3f7;
  border: 1px solid #d4dde5;
  border-radius: 50%;
  color: #314554;
  display: inline-flex;
  font-size: 12px;
  font-weight: 700;
  height: 28px;
  justify-content: center;
  width: 28px;
}

.ordered-sequence-identity {
  display: flex;
  flex-direction: column;
  min-width: 0;
}

.ordered-sequence-identity strong {
  color: #1f2933;
  font-size: 14px;
  line-height: 1.25;
  overflow-wrap: anywhere;
}

.ordered-sequence-identity span {
  color: #61707c;
  font-size: 12px;
  font-weight: 600;
  line-height: 1.3;
}

.ordered-sequence-actions {
  display: flex;
  gap: 4px;
}

.ordered-sequence-action {
  align-items: center;
  background: #ffffff;
  border: 1px solid #c8d0d7;
  border-radius: 4px;
  color: #314554;
  display: inline-flex;
  height: 30px;
  justify-content: center;
  padding: 0;
  width: 30px;
}

.ordered-sequence-action:hover:not(:disabled),
.ordered-sequence-action:focus:not(:disabled) {
  border-color: #226597;
  color: #226597;
  outline: 0;
}

.ordered-sequence-action:disabled {
  color: #a8b0b7;
  cursor: default;
  opacity: 0.55;
}

.ordered-sequence-remove:hover:not(:disabled),
.ordered-sequence-remove:focus:not(:disabled) {
  border-color: #a83a3a;
  color: #a83a3a;
}

.ordered-sequence-source {
  border-color: #b7c0c9;
  border-radius: 6px;
  font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", monospace;
  font-size: 12px;
  line-height: 1.45;
  min-height: 58px;
  resize: vertical;
  white-space: pre;
}

.data-entry-script-preview {
  background: #17212b;
  border-radius: 6px;
  color: #f3f7fa;
  font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", monospace;
  font-size: 12px;
  line-height: 1.45;
  margin: 0;
  max-height: 450px;
  min-height: 360px;
  overflow: auto;
  overflow-wrap: normal;
  padding: 14px;
  white-space: pre;
}

.validation-copy {
  height: 1px;
  left: -10000px;
  opacity: 0;
  overflow: hidden;
  pointer-events: none;
  position: absolute;
  top: auto;
  width: 1px;
}

.run-actions-bar {
  display: flex;
  gap: 12px;
  justify-content: flex-end;
  margin-top: 80px;
  max-width: 1017px;
}

.run-action-button {
  min-width: 120px;
}

@media (max-width: 920px) {
  .new-run-content {
    padding: 16px 18px 70px;
  }

  .new-run-heading {
    align-items: flex-start;
    display: grid;
    gap: 12px;
  }

  .run-type-selection-grid {
    grid-template-columns: 1fr;
  }

  .workflow-selection-grid {
    grid-template-columns: 1fr;
  }

  .run-type-switcher-row,
  .workflow-stepper,
  .workflow-step {
    align-items: flex-start;
    flex-direction: column;
  }

  .run-type-switcher-row {
    grid-template-columns: 1fr;
  }

  .workflow-step + .workflow-step::before {
    display: none;
  }

  .manual-field-row,
  .output-record-row,
  .ordered-sequence-toolbar {
    grid-template-columns: 1fr;
  }

  .ordered-sequence-add {
    justify-self: start;
  }

  .output-record-extension {
    text-align: left;
  }
}
</style>
