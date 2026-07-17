import re
from pathlib import Path


CREATE_RUN = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "components"
    / "Pages"
    / "CreateRun.vue"
)
HOME = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "components"
    / "Pages"
    / "Home.vue"
)
WORKFLOW_RUN = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "components"
    / "Pages"
    / "WorkflowRun.vue"
)
RESOURCE_SETTINGS = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "components"
    / "blocks"
    / "RunResourceSettings.vue"
)
VIEW_RUN = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "components"
    / "Pages"
    / "ViewRun.vue"
)
ROUTER = Path(__file__).resolve().parents[1] / "src" / "router.js"
USER_STORAGE = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "components"
    / "overlay"
    / "UserStorage.vue"
)
EPOLYSCAT_SERVICE = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "service"
    / "epolyscat-service.js"
)
SETTINGS_STORE = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "store"
    / "modules"
    / "settings.store.js"
)
INPUT_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "utils"
    / "epolyscat-input-script.js"
)


def _source():
    return CREATE_RUN.read_text()


def _home_source():
    return HOME.read_text()


def _workflow_source():
    return WORKFLOW_RUN.read_text()


def _resource_settings_source():
    if not RESOURCE_SETTINGS.exists():
        return ""
    return RESOURCE_SETTINGS.read_text()


def _input_script_source():
    return INPUT_SCRIPT.read_text()


def _view_run_source():
    return VIEW_RUN.read_text()


def _router_source():
    return ROUTER.read_text()


def _user_storage_source():
    return USER_STORAGE.read_text()


def _epolyscat_service_source():
    return EPOLYSCAT_SERVICE.read_text()


def _route_block(path):
    source = _router_source()
    pattern = re.compile(
        r"\{\s*path:\s*['\"]" + re.escape(path) + r"['\"],(?P<body>.*?)\n\s*\}",
        re.DOTALL,
    )
    match = pattern.search(source)
    assert match, f"Route {path} is missing"
    return match.group("body")


def test_create_run_page_uses_new_run_layout_hooks():
    source = _source()
    resource_source = _resource_settings_source()

    assert 'class="new-run-page"' in source
    assert 'class="new-run-content"' in source
    assert 'class="input-files-panel"' in source
    assert 'class="data-entry-card"' in source
    assert "<RunResourceSettings" in source
    assert 'class="resource-settings-grid"' in resource_source


def test_create_run_selects_run_type_application_and_required_files_inline():
    source = _source()

    expected_hooks = [
        'class="run-type-selection-grid"',
        'class="run-selection-column run-type-column"',
        "Run Type",
        "{{ activeRunTypeTitle }}",
        "Required Files",
        "selectedRunType",
        "runTypeOptions",
        "activeRunTypeApplications",
        "activeRequiredFiles",
        'v-for="runType in runTypeOptions"',
        'v-for="application in activeRunTypeApplications"',
        'v-for="file in activeRequiredFiles"',
        "selectRunType(runType.id)",
        "selectRunApplication(application.id)",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert 'v-if="selectedRunType !== \'workflow\'"' in source


def test_create_run_replicates_legacy_descriptions_panel():
    source = _source()

    expected_hooks = [
        "Descriptions",
        "showDescriptions",
        "new-run-descriptions",
        "new-run-description-list",
        "newRunDescriptions()",
        "descriptionPath()",
        "descriptions[item]",
        'import { descriptions } from "@/fileData";',
        "this.$store.commit(\"settings/setPreference\"",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_create_run_hides_file_sections_when_selection_has_no_required_files():
    source = _source()

    expected_hooks = [
        "hasRequiredFiles()",
        '<section v-if="hasRequiredFiles" class="input-files-panel"',
        '<section v-if="hasRequiredFiles && isEPolyScatScriptInput" class="data-entry-section"',
        "inpcContent: !this.hasRequiredFiles ? true",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_create_run_catalog_maps_module_utility_workflow_to_required_files():
    source = _source()

    expected_catalog_tokens = [
        "Modules",
        "Utilities",
        "Workflows",
        "Gaussian16",
        "OpenMolcas",
        "ePolyScat",
        "Data_Gen",
        "ePolyScat_Run",
        "Analysis",
        "CnvMath",
        "MoldenMerge",
        "Gaussian_Input",
        "Molcas_Input",
        "ePolyScat_Input_Data",
        "ePolyscat_Input_File",
        "molden.dat",
    ]

    for token in expected_catalog_tokens:
        assert token in source


def test_create_run_uses_legacy_utility_names_with_full_descriptions():
    source = _source()

    expected_utility_labels = [
        'label: "CnvMath"',
        'description: "ConvertToMathematica"',
        'label: "CnvMatLab"',
        'description: "ConvertToMatlab"',
        'label: "CnvLinFull"',
        'description: "Compute Ph.ion.diff.Xsec."',
        'label: "MoldenMerge"',
        'description: "Merge Molden Data Files"',
        'label: "NRFPAD"',
        'description: "N photon RFPAD"',
        'label: "Cube2igor"',
        'description: "G16CubeToIGOR Plots"',
    ]

    for token in expected_utility_labels:
        assert token in source


def test_analysis_utilities_declare_manual_specific_required_inputs():
    source = _source()
    service_source = _epolyscat_service_source()

    expected_create_run_tokens = [
        'id: "CnvLinFull"',
        'name: "DumpOut"',
        'description: "DumpIdy output file"',
        'id: "CnvMath"',
        'name: "BendOrient_Output"',
        'id: "CnvMatLab"',
        'id: "NRFPAD"',
        'name: "Cross_Section_Input_File"',
        'id: "Cube2igor"',
        'name: "Cube_Output"',
    ]
    expected_service_tokens = [
        '"name": "DumpOut"',
        '"value": "CnvLinFull"',
        '"name": "BendOrient_Output"',
        '"value": "CnvMath"',
        '"value": "CnvMatLab"',
        '"name": "Cross_Section_Input_File"',
        '"value": "NRFPAD"',
        '"name": "Cube_Output"',
        '"value": "Cube2igor"',
    ]

    for token in expected_create_run_tokens:
        assert token in source
    for token in expected_service_tokens:
        assert token in service_source


def test_create_run_required_files_are_status_links_to_input_files():
    source = _source()

    expected_hooks = [
        'ref="inputFilesPanel"',
        "required-file-status",
        "required-file-status-ready",
        "required-file-status-missing",
        "uploadedRequiredFiles(file.name).length > 0",
        "!file.generatedByNewRun",
        'v-on:click="selectRequiredFile(file.name)"',
        "selectRequiredFile(fileName)",
        "scrollToInputFiles()",
        "this.$refs.inputFilesPanel.scrollIntoView",
    ]

    for hook in expected_hooks:
        assert hook in source

    forbidden_hooks = [
        "required-file-actions",
        "new-run-required-file-storage",
        "openRequiredFileStorage",
        "onRequiredStorageFilesSelected",
        "onRequiredComputerFilesSelected",
        "addRequiredFilesToInput",
    ]

    for hook in forbidden_hooks:
        assert hook not in source


def test_create_run_required_file_status_ignores_generated_placeholders():
    source = _source()
    method = re.search(
        r"uploadedRequiredFiles\(fileName\) \{(?P<body>.*?)\n    \},",
        source,
        re.DOTALL,
    )

    assert method, "uploadedRequiredFiles method is missing"
    assert "!file.deleted" in method.group("body")
    assert "!file.generatedByNewRun" in method.group("body")


def test_create_run_workflow_application_picker_is_separated_from_stage_options():
    source = _source()

    expected_hooks = [
        "workflow-stepper",
        "workflow-step",
        "workflow-step-index",
        "workflow-step-label",
        "workflow-stage-application-panel",
        "workflow-stage-application-options",
        "workflow-stage-application-option",
        "activeWorkflowApplicationTitle",
        "selectedWorkflowApplicationId",
        "selectWorkflowApplication(application.id)",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "workflow-application-picker" not in source


def test_create_run_input_files_section_owns_upload_actions():
    source = _source()

    expected_hooks = [
        "input-file-actions",
        "Select file from storage",
        "Select file from computer",
        "<UserStorage",
        'id="new-run-input-file-storage"',
        ':canSelectMultiple="activeInputAllowsMultiple"',
        '@filesSelected="onInputStorageFilesSelected"',
        'type="file"',
        'v-on:change="onInputComputerFilesSelected"',
        "addFilesToActiveInput",
        'this.$store.commit("input/ADD_TO_INPUT_FILE"',
        "input-file-remove",
        'v-on:click="removeInputFile(file)"',
        "removeInputFile(file)",
        'this.$store.commit("input/REMOVE_FILE"',
        "filename: file.name",
        "inputFileName: this.activeInputFile.inputFileName",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_create_run_save_payload_uses_inline_run_selection():
    source = _source()

    expected_payload_tokens = [
        "selectedRunMode()",
        "selectedModuleApplication()",
        "selectedWorkflowStage()",
        "selectedWorkflowApplication()",
        "selectedUtilityApplication()",
        "runMode: this.selectedRunMode",
        "moduleApplication: this.selectedModuleApplication",
        "workflowStage: this.selectedWorkflowStage",
        "workflowApplication: this.selectedWorkflowApplication",
        "utilityApplication: this.selectedUtilityApplication",
        "isWorkflowPlan: this.selectedRunType === \"workflow\"",
        "dataGenerationApplication: this.workflowStageSelections.Data_Gen",
        "analysisApplications: [this.workflowStageSelections.Analysis",
        "plannedStageIds: [\"Data_Gen\", \"ePolyScat_Run\", \"Analysis\", \"Visualization\"]",
    ]

    for token in expected_payload_tokens:
        assert token in source


def test_create_run_workflow_submit_guards_later_steps_until_prior_steps_complete():
    source = _source()
    method = re.search(
        r"validateWorkflowSubmitStep\(\) \{(?P<body>.*?)\n    \},",
        source,
        re.DOTALL,
    )

    expected_hooks = [
        "validateWorkflowSubmitStep()",
        "firstWorkflowSubmitStage()",
        "this.selectedApplicationId === firstSubmitStage.id",
        "this.selectRunApplication(firstSubmitStage.id)",
        'ref="newRunTopAnchor"',
        'id="new-run-top-anchor"',
        '<b-modal',
        'v-model="showWorkflowStepGuardModal"',
        'return-focus="#new-run-top-anchor"',
        'v-on:hidden="onWorkflowStepGuardModalHidden"',
        "workflowStepGuardMessage",
        "showWorkflowStepGuardModal: false",
        "workflowStepGuardPendingScroll: false",
        "this.showWorkflowStepGuardModal = true",
        "this.workflowStepGuardPendingScroll = true",
        "onWorkflowStepGuardModalHidden()",
        "window.requestAnimationFrame",
        "this.scrollToNewRunTop()",
        "scrollToNewRunTop()",
        'ref="newRunContent"',
        "document.activeElement.blur",
        'closest(".overflow-auto")',
        "scrollContainer.scrollTo",
        "top: 0",
        "this.$refs.newRunContent.scrollIntoView",
        "Please complete",
        "before submitting this workflow.",
        "if (!this.validateWorkflowSubmitStep()) return;",
    ]

    for token in expected_hooks:
        assert token in source

    assert method, "validateWorkflowSubmitStep method is missing"
    assert 'eventBus.$emit("error"' not in method.group("body")
    assert "this.$bvModal.msgBoxOk" not in method.group("body")


def test_create_run_data_entry_models_epolyscat_input_script_sections():
    create_source = _source()
    source = create_source + _input_script_source()

    expected_tokens = [
        "input_data",
        "input_file (.inp)",
        "Grid / Expansion",
        "State Definitions",
        "Potentials",
        "Energies / Partial Waves",
        "Outputs",
        "dataEntrySections",
        "activeDataEntrySection",
        "data-entry-script-preview",
        "buildEPolyScatInputScript",
        "LMax",
        "EMax",
        "ScatEng",
        "FileName",
    ]

    for token in expected_tokens:
        assert token in source

    assert "placeholder-tab-panel" not in create_source
    assert "This tab is waiting for the next screenshot." not in create_source
    assert "input_data (.inp)" not in create_source
    assert "Command Sequence" not in create_source
    assert "command-sequence" not in create_source


def test_create_run_embeds_workflow_stepper_when_workflows_selected():
    source = _source()

    expected_hooks = [
        'v-else',
        'class="run-type-switcher-row"',
        'class="workflow-stepper"',
        'class="workflow-step"',
        'class="workflow-step-index"',
        'v-for="(application, index) in activeRunTypeApplications"',
        "{{ index + 1 }}",
        "Data_Gen",
        "ePolyScat_Run",
        "Analysis",
        "Post-processing",
        "Visualization",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "new-run-workflow-stepper" not in source
    assert "workflowStages" not in source
    assert "selectWorkflowStage(stage)" not in source
    assert "wizard" not in source.lower()


def test_home_create_run_button_opens_new_run_catalog():
    source = _home_source()

    expected_hooks = [
        "goToNewRun",
        '@click="goToNewRun"',
        'this.$router.push("/runs/new")',
    ]

    for hook in expected_hooks:
        assert hook in source

    assert re.search(
        r"<b-button[^>]*@click=\"goToNewRun\"[^>]*>\s*Create new run\s*</b-button>",
        source,
        re.DOTALL,
    )

    assert "create-run-type-modal" not in source
    assert "chooseWorkflowRun" not in source


def test_create_run_page_does_not_auto_show_entry_choice_modal():
    source = _source()

    forbidden_hooks = [
        "create-run-type-modal",
        "showRunTypeModal",
        "chooseModuleRun",
        "chooseWorkflowRun",
        "skipCreateRunTypeModal",
    ]

    for hook in forbidden_hooks:
        assert hook not in source


def test_create_run_and_workflow_pages_offer_run_type_switches():
    create_source = _source()
    workflow_source = _workflow_source()

    assert "run-type-switcher-row" in create_source
    assert "run-selection-column run-type-column" in create_source
    assert "selectRunType(runType.id)" in create_source
    assert "Switch to Workflow Run" not in create_source
    assert "switchToModuleRun" in workflow_source
    assert "Switch to Module Run" in workflow_source


def test_create_run_output_section_tracks_manual_file_outputs():
    source = _source() + _input_script_source()

    expected_output_tokens = [
        "outputDefinitions",
        "MatrixElements",
        "PlotData",
        "DumpOut",
        "OrientData",
        "ViewOrb",
        "ViewOrbGeom",
        ".idy",
        ".dat",
        "fileType",
        "fileName",
    ]

    for token in expected_output_tokens:
        assert token in source


def test_create_run_table_view_does_not_duplicate_file_command_sequence():
    source = _source()

    assert "ePolyScatInputScript" in source
    assert "buildEPolyScatInputScript" in source
    assert "commandSequence" not in source
    assert "command-sequence-list" not in source
    assert "command-step" not in source
    assert 'type: "commands"' not in source
    assert "targetStateInputColumns" not in source


def test_create_run_omits_save_to_dropdown_from_resource_settings():
    source = _source()

    removed_hooks = [
        ':show-save-target="true"',
        ':save-target="saveTarget"',
        ':save-target-options="saveTargetOptions"',
        'v-on:saveTargetSelected="selectSaveTarget"',
        "saveTarget:",
        "saveTargetOptions",
        "selectSaveTarget(option)",
        "[Save To]",
    ]

    for hook in removed_hooks:
        assert hook not in source


def test_create_run_preserves_resource_and_submission_controls():
    source = _source()
    resource_source = _resource_settings_source()

    assert "GroupResourceProfileService.list" in resource_source
    assert "ComputeResourceService.names" in resource_source
    assert "ApplicationDeploymentService.list" in resource_source
    assert "adpf-queue-settings-editor" in resource_source
    assert 'v-on:click="onSave(false)"' in source


def test_create_run_filters_resources_for_the_selected_application_module():
    create_source = _source()
    resource_source = _resource_settings_source()
    settings_source = SETTINGS_STORE.read_text()

    assert ':application-module-id="resourceApplicationModuleId"' in create_source
    assert "resourceApplicationModuleId()" in create_source
    assert "gaussian16ApplicationModuleId" in create_source
    assert "openmolcasApplicationModuleId" in create_source
    assert "applicationModuleId" in resource_source
    assert "this.applicationModuleId" in resource_source
    assert "GAUSSIAN16_APPLICATION_ID" in settings_source
    assert "OPENMOLCAS_APPLICATION_ID" in settings_source
    assert "gaussian16ApplicationModuleId" in settings_source
    assert "openmolcasApplicationModuleId" in settings_source


def test_create_run_and_workflow_reuse_same_resource_settings_component():
    create_source = _source()
    workflow_source = _workflow_source()

    assert "RunResourceSettings" in _resource_settings_source()
    assert "RunResourceSettings" in create_source
    assert "RunResourceSettings" in workflow_source
    assert "@/components/blocks/RunResourceSettings" in create_source
    assert "@/components/blocks/RunResourceSettings" in workflow_source
    assert "@/components/blocks/RunResource\"" not in workflow_source
    assert "<RunResource " not in workflow_source


def test_workflow_page_owns_stepper_and_application_selection():
    source = _workflow_source()

    expected_hooks = [
        'class="workflow-run-page"',
        "workflow-stepper",
        "Workflow/{{ activeStageLabel }}",
        "Data Generation",
        "ePolyScat Run",
        "Post-processing",
        "Visualization",
        "Gaussian16",
        "OpenMolcas",
        "Select file from storage",
        "Select file from computer",
        "workflowApplication",
        "runMode: \"workflow\"",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "input_data (.inp)" not in source


def test_view_run_page_renders_module_and_workflow_read_models():
    source = _view_run_source()

    expected_hooks = [
        'class="view-run-page"',
        "presentation.mode",
        "Modules/EPOLYSCAT_DMAT",
        "Workflow/STGF",
        "file-selector",
        "target-states-matrix",
        "run-code-viewer",
        "plot-panel",
        "RunResource",
        "fetchRunPresentation",
        "workflow-status-panel",
        "workflowStatusLabel",
        "activeWorkflowStage",
        "Current step",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_view_run_workflow_panel_distinguishes_child_step_from_parent_workflow():
    source = _view_run_source()

    expected_hooks = [
        'v-for="item in workflowPanelItems"',
        "{{ item.label }}",
        "{{ item.value }}",
        "isWorkflowChildRun",
        "Step status",
        "Workflow position",
        "Next step",
        "stepStatusLabel",
        "workflowPositionLabel",
        "nextWorkflowStage",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert '<span class="workflow-status-label">Workflow status</span>' not in source
    assert '<span class="workflow-status-label">Current step</span>' not in source
    assert "workflowPanelLink" not in source
    assert "workflow-status-link" not in source


def test_view_run_fetches_real_file_catalogs_and_selected_file_content():
    source = _view_run_source()

    expected_hooks = [
        "InputService.fetchOutputs",
        "PlotService.getViewables",
        "PlotService.getInputFiles",
        "RunService.fetchViewableContent",
        "fetchSelectedFileContent",
        "selectedFileContent",
        "filePreviewLoading",
        "filePreviewError",
        "watch:",
        "selectedFile()",
        "resolvedFileGroups",
        "encodeURIComponent(filename)",
    ]

    for hook in expected_hooks:
        assert hook in source or hook in (
            Path(__file__).resolve().parents[1]
            / "src"
            / "service"
            / "epolyscat-service.js"
        ).read_text()


def test_view_run_parent_workflow_uses_active_child_for_file_catalog():
    source = _view_run_source()

    expected_hooks = [
        "fileCatalogRunId",
        'this.presentation.mode === "workflow"',
        "!this.isWorkflowChildRun",
        "this.activeWorkflowStage.child_run_id",
        "PlotService.getViewables({ runId: this.fileCatalogRunId })",
        "PlotService.getInputFiles({ runId: this.fileCatalogRunId })",
        "InputService.fetchOutputs(this.fileCatalogRunId)",
        "runId: this.fileCatalogRunId",
        "runIds: [this.fileCatalogRunId]",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_view_run_pending_workflow_step_links_to_configure_page():
    source = _view_run_source()

    expected_hooks = [
        "workflowStageRoute(stage, index)",
        "workflowStageIsConfigurable(stage, index)",
        "workflow-step-configure-link",
        "workflowConfigureStepLink(stage, index)",
        "workflowChildRunId=${stage.child_run_id}",
        "workflowParentRunId=${this.run.id}",
        "withOutputsFrom=${previousStage.child_run_id}",
        'previousStage.state === "complete"',
        'previousStage.status === "complete"',
        "return `/runs/${stage.child_run_id}`",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_view_run_can_continue_an_eligible_completed_run_into_workflow():
    source = _view_run_source()
    service_source = _epolyscat_service_source()

    expected_view_hooks = [
        "Continue in Workflow",
        "workflowContinuation",
        "workflowContinuationLoading",
        "loadWorkflowContinuation",
        "continueInWorkflow",
        "RunService.fetchWorkflowContinuation",
        "RunService.continueWorkflow",
        "workflowChildRunId=${continuation.nextChildRunId}",
        "workflowParentRunId=${continuation.workflowParentRunId}",
        "withOutputsFrom=${continuation.sourceRunId}",
        "workflowStageStatusText(stage)",
        'previousStage.status === "imported"',
        'previousStage.status === "not_included"',
    ]
    expected_service_hooks = [
        "fetchWorkflowContinuation",
        "continueWorkflow",
        "workflow_continuation/",
        "workflow_parent_run_id",
        "next_child_run_id",
        "source_run_id",
    ]

    for hook in expected_view_hooks:
        assert hook in source
    for hook in expected_service_hooks:
        assert hook in service_source


def test_create_run_can_edit_existing_workflow_child_step():
    source = _source()

    expected_hooks = [
        "Workflow Step",
        "isWorkflowChildEdit",
        "workflowChildRunId()",
        "workflowParentRunId()",
        "this.$route.query.workflowChildRunId",
        "applyWorkflowChildRun(sourceRun)",
        'this.selectedRunType = "workflow"',
        "this.selectedApplicationId = stageId",
        'this.$store.dispatch("run/loadInputs", { runId: workflowChildRunId })',
        "const outputsRunId = parseInt(this.outputsRunId)",
        "const outputSourceRun = await this.$store.dispatch(\"run/fetchRun\", { runId: outputsRunId })",
        "this.loadWorkflowPreviousOutputsIntoActiveInputs(outputsRunId, outputSourceRun)",
        'this.$store.commit("input/SET_PATH", { path: this.selectedRunPath })',
        "addWorkflowOutputToInput(inputFileName, outputFile)",
        "loadWorkflowPreviousOutputsIntoActiveInputs(outputsRunId, sourceRun = null)",
        "InputService.fetchWorkflowOutputBinding({",
        "targetStageId: this.selectedApplicationId",
        "targetApplicationId: this.activeWorkflowApplicationId",
        "requiredFileName",
        'if (backendBinding.status === "ready")',
        "backendBinding.outputFile",
        "InputService.fetchOutputs(outputsRunId)",
        "buildWorkflowOutputInputBinding",
        "workflowApplicationIdFromRun(sourceRun)",
        "applyWorkflowDataEntryValues(binding.dataEntryValues)",
        "this.selectedApplicationId === \"ePolyScat_Run\"",
        "targetApplicationId: this.activeWorkflowApplicationId",
        "if (this.isWorkflowChildEdit) {",
        'this.$store.dispatch("run/updateRun", runPayload)',
        "parent_run_id: this.run.parentRunId || this.workflowParentRunId",
        'return this.isWorkflowChildEdit ? this.selectedApplicationId : ""',
        "router.push(`/runs/${this.run.parentRunId || this.workflowParentRunId || this.run.id}`)",
    ]

    for hook in expected_hooks:
        assert hook in source

    service_source = _epolyscat_service_source()
    expected_service_hooks = [
        "fetchWorkflowOutputBinding",
        "workflow_output_binding/",
        "targetStageId",
        "targetApplicationId",
        "requiredFileName",
        "inputFileName: data.input_file_name",
        "const outputFile = data.selected",
        "dataEntryValues: data.data_entry_values",
    ]
    for hook in expected_service_hooks:
        assert hook in service_source


def test_view_run_plot_panel_uses_plot_service_not_static_svg():
    source = _view_run_source()

    expected_hooks = [
        "PlotService.plotSelectedRuns",
        "plotImageUrl",
        "plotFileOptions",
        "canCreatePlot",
        "createPlot",
        'v-on:click="createPlot"',
    ]

    for hook in expected_hooks:
        assert hook in source

    assert '<svg viewBox="0 0 260 160"' not in source
    assert "<polyline" not in source


def test_view_run_plot_panel_is_gated_by_run_type_and_plottable_outputs():
    source = _view_run_source()

    expected_hooks = [
        'v-if="showPlotPanel"',
        ":class=\"{ 'with-plot-panel': showPlotPanel }\"",
        "showPlotPanel",
        "isPlotCapableRun",
        "activeRunApplicationForPlot",
        "plottableOutputFiles",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert '<aside class="plot-panel" v-if="presentation.mode === \'workflow\'">' not in source
    assert ':class="{ \'workflow-layout\': presentation.mode === \'workflow\' }"' not in source


def test_view_run_plot_dropdown_only_uses_plottable_files():
    source = _view_run_source()

    expected_hooks = [
        "plottableFileNames",
        "presentation.plottable_file_names",
        "isPlottableFile(file)",
        "plottableOutputFileForName",
        "const files = this.plottableOutputFiles",
        "const outputFile = this.plottableOutputFileForName(this.plotForm.file)",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "const files = this.outputFiles.filter(file => this.fileDataProductURI(file));" not in source
    assert "const outputFile = this.outputFileForName(this.plotForm.file);" not in source


def test_view_run_file_selector_uses_resolved_real_file_groups():
    source = _view_run_source()

    assert 'v-for="group in resolvedFileGroups"' in source
    assert "presentationGroupFiles" in source
    assert "outputFileForName" in source


def test_view_run_reads_selected_input_file_from_run_inputs():
    source = _view_run_source()

    expected_hooks = [
        "runInputFiles",
        "inputFileForName",
        "const inputFile = this.inputFileForName(this.selectedFile)",
        "InputService.fetchFileContents(inputFile)",
        "this.run.inputs",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_view_run_header_deduplicates_status_badges():
    source = _view_run_source()

    expected_hooks = [
        "statusBadges",
        'v-for="badge in statusBadges"',
        "badge.label",
        "badge.value",
        "primaryStatusLabel",
        'return "Step";',
        'return "Workflow";',
        "jobStatus !== runStatus",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "{{ run.displayStatus || run.status }}" not in source
    assert "{{ run.jobStatus }}" not in source
    assert 'badges.push({ label: "Run", value: runStatus })' not in source


def test_view_run_target_matrix_uses_resolved_input_and_output_files():
    source = _view_run_source()

    expected_hooks = [
        "targetInputColumns",
        "targetOutputFiles",
        'v-for="(column, columnIndex) in targetInputColumns"',
        'v-for="file in targetOutputFiles"',
        "chunkFileNames",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert 'v-for="(column, columnIndex) in presentation.target_states.input_columns"' not in source
    assert 'v-for="file in presentation.target_states.output_files"' not in source


def test_view_run_target_output_uses_same_real_files_as_file_selector():
    source = _view_run_source()

    expected_hooks = [
        "resolvedOutputFiles",
        "return this.resolvedOutputFiles;",
        "const outputNames = this.resolvedOutputFiles;",
        "...viewableNames.filter(filename => inputNames.indexOf(filename) < 0)",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "const actualOutputFiles = this.uniqueFileNames(this.outputFiles.map(file => this.fileDisplayName(file)));" not in source
    assert '...this.presentationGroupFiles("Outputs")' not in source


def test_view_run_resolved_output_files_excludes_input_names():
    source = _view_run_source()

    expected_hooks = [
        "const outputFileNames = this.outputFiles.map(file => this.fileDisplayName(file));",
        "...outputFileNames.filter(filename => inputNames.indexOf(filename) < 0)",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_view_run_target_matrix_uses_compact_adaptive_layout():
    source = _view_run_source()

    expected_hooks = [
        "align-items: start;",
        "repeat(auto-fit, minmax(260px, 1fr))",
        "repeat(auto-fit, minmax(150px, 1fr))",
        'this.presentationGroupFiles("Inputs")',
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "grid-template-columns: minmax(0, 506px) minmax(0, 193px);" not in source


def test_view_run_code_viewer_is_file_view_only():
    source = _view_run_source()
    toolbar = re.search(
        r"<div class=\"code-viewer-toolbar\">(?P<body>.*?)</div>",
        source,
        re.DOTALL,
    )

    assert toolbar, "code viewer toolbar is missing"
    assert "<b-icon" not in toolbar.group("body")
    assert 'icon="list-ul"' not in source
    assert 'icon="card-list"' not in source


def test_view_run_file_selector_lives_in_file_view_toolbar():
    source = _view_run_source()
    code_viewer = re.search(
        r"<section class=\"run-code-viewer\">(?P<body>.*?)</section>",
        source,
        re.DOTALL,
    )

    assert code_viewer, "run-code-viewer section is missing"
    assert 'class="file-control-row"' in code_viewer.group("body")
    assert 'for="run-file-selector"' in code_viewer.group("body")
    assert 'id="run-file-selector"' in code_viewer.group("body")
    assert re.search(
        r"<section class=\"view-run-main\">(?P<body>.*?)<section class=\"target-states-matrix\"",
        source,
        re.DOTALL,
    )
    assert 'class="file-control-row"' not in re.search(
        r"<section class=\"view-run-main\">(?P<body>.*?)<section class=\"target-states-matrix\"",
        source,
        re.DOTALL,
    ).group("body")


def test_view_run_main_content_is_centered_when_plot_panel_is_hidden():
    source = _view_run_source()

    expected_hooks = [
        ":class=\"{ 'with-plot-panel': showPlotPanel }\"",
        ".view-run-layout:not(.with-plot-panel) .view-run-main",
        "justify-self: center;",
        "max-width: 760px;",
        "width: 100%;",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_router_splits_create_workflow_and_view_run_pages():
    source = _router_source()

    assert "WorkflowRun" in source
    assert "ViewRun" in source
    assert "path: '/runs/new/workflow'" in source
    assert "name: 'WorkflowRun'" in source
    assert "name: 'Run'" in source


def test_create_run_reuses_existing_input_store_binding():
    source = _source()

    expected_hooks = [
        'this.$store.dispatch("input/fetchPathLabels")',
        'this.$store.commit("input/SET_PATH"',
        'this.$store.commit("input/ADD_TO_INPUT_FILE"',
        "selectedRunPath",
        "generatedInputFileName",
        "syncGeneratedInputFile",
        "activeInputFile",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert 'inputFileName: "ePolyScat_Input_Data"' in source
    assert 'inputFileName: "ePolyscat_Input_File"' in source
    assert 'generatedFileName: "input_data"' in source
    assert 'generatedFileName: "input_file.inp"' in source


def test_input_file_tabs_sync_all_required_old_input_files_before_save():
    source = _source()

    assert "syncAllGeneratedInputFiles" in source
    assert "this.inputFiles.forEach" in source
    assert "this.syncAllGeneratedInputFiles();" in source


def test_workflow_file_inputs_expose_upload_and_gateway_storage_actions():
    source = _workflow_source()

    expected_hooks = [
        "workflow-file-actions",
        "Select file from storage",
        "Select file from computer",
        "<UserStorage",
        '@filesSelected="onStorageFilesSelected"',
        'type="file"',
        'v-on:change="onComputerFileSelected"',
        "selectedWorkflowFile",
        "selectedWorkflowFileName",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_workflow_upload_uses_single_input_file_selector():
    source = _workflow_source()

    assert 'type="file"' in source
    assert "\n                multiple\n" not in source


def test_workflow_storage_file_picker_is_single_file_selection():
    source = _workflow_source()
    user_storage_source = _user_storage_source()

    assert ':canSelectMultiple="false"' in source
    assert 'modalTitle="Select workflow input file"' in source
    assert "canSelectMultiple" in user_storage_source
    assert ':canSelectMultiple="allowMultipleSelection"' in user_storage_source
    assert "allowMultipleSelection()" in user_storage_source
    assert 'return this.selectionType !== "folder";' in user_storage_source
    assert "modalTitle" in user_storage_source


def test_workflow_run_uses_canonical_airavata_workflow_inputs():
    source = _workflow_source()

    expected_tokens = [
        'activeStageId: "Data_Gen"',
        '{ id: "Data_Gen", label: "Data Generation"',
        '{ id: "ePolyScat_Run", label: "ePolyScat Run"',
        '{ id: "Analysis", label: "Post-processing"',
        '{ id: "Visualization", label: "Visualization"',
        "activeWorkflowInputName()",
        "`Application_Workflow = ${this.activeStageId}`",
        "`Data_Gen = ${this.workflowApplication}`",
        '{ type: "radio input", name: "Application_Workflow", value: this.activeStageId }',
        '{ type: "radio input", name: "Data_Gen", value: this.workflowApplication }',
        "name: this.activeWorkflowInputName",
    ]

    for token in expected_tokens:
        assert token in source

    assert "Workflow_Application" not in source


def test_new_run_workflow_splits_post_processing_from_local_visualization():
    source = _source()

    assert 'id: "Analysis"' in source
    assert 'label: "Post-processing"' in source
    assert 'id: "Visualization"' in source
    assert 'label: "Visualization"' in source
    assert 'localOnly: true' in source
    assert 'plannedStageIds: ["Data_Gen", "ePolyScat_Run", "Analysis", "Visualization"]' in source


def test_view_run_opens_ready_visualization_on_post_processing_outputs():
    source = _view_run_source()

    assert "stage.local_only" in source
    assert "stage.source_child_run_id" in source
    assert 'query: { visualize: "1" }' in source
    assert 'stage.status === "ready"' in source
    assert 'this.$route.query.visualize === "1"' in source
    assert "this.visualizationMode" in source
    assert 'return "Workflow/Visualization"' in source


def test_workflow_file_source_handlers_keep_computer_and_storage_paths_separate():
    source = _workflow_source()
    computer_method = re.search(
        r"onComputerFileSelected\(event\) \{(?P<body>.*?)\n    \},",
        source,
        re.DOTALL,
    )
    storage_method = re.search(
        r"onStorageFilesSelected\(files\) \{(?P<body>.*?)\n    \},",
        source,
        re.DOTALL,
    )
    assert computer_method, "onComputerFileSelected method is missing"
    assert storage_method, "onStorageFilesSelected method is missing"

    assert "file.isFromComputer = true" in computer_method.group("body")
    assert "isFromComputer: false" in storage_method.group("body")


def test_input_file_tabs_drive_visible_data_entry_content():
    source = _source()

    assert "isEPolyScatScriptInput()" in source
    assert 'activeInputFile.id === "input_file"' in source
    assert '<section v-if="hasRequiredFiles && isEPolyScatScriptInput" class="data-entry-section">' in source


def test_data_entry_view_mode_changes_visible_layout():
    source = _source()

    expected_hooks = [
        ':class="dataEntryCardClasses"',
        "dataEntryCardClasses()",
        'dataEntryViewMode: "table"',
        "dataEntryViewMode === 'file'",
        "dataEntryViewMode === 'table'",
        'v-on:click="selectDataEntryViewMode(\'file\')"',
        'v-on:click="selectDataEntryViewMode(\'table\')"',
        'title="File view"',
        'title="Table view"',
        'v-if="dataEntryViewMode === \'table\'"',
        'v-else class="data-entry-file-view"',
        ".data-entry-card.table-mode",
        ".data-entry-file-view",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "data-entry-script-panel" not in source
    assert "data-entry-editor-grid-list" not in source
    assert ".data-entry-card.list-mode" not in source


def test_data_entry_file_view_is_editable_and_syncs_to_table_values():
    source = _source()

    expected_hooks = [
        "<b-form-textarea",
        'v-model="inpcContent"',
        'v-on:input="onFileViewInput"',
        'v-on:input="onTableFieldInput(field.key, $event)"',
        'v-on:input="onTableFieldInput(output.valueKey, $event)"',
        "selectDataEntryViewMode(mode)",
        "applyInpcContentToDataEntryValues",
        "releaseFileContentApplyLock",
        "parseEPolyScatInputScript",
        "isApplyingFileContent",
        "patchInpcContentFromTableValues",
        "patchEPolyScatInputScript",
        "syncInpcContentToActiveInputFile",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_uploaded_file_table_uses_blank_missing_values_with_default_placeholders():
    source = _source()

    expected_hooks = [
        "dataEntryRecommendedValues",
        "blankDataEntryValues()",
        "dataEntryPlaceholder(field.key)",
        "dataEntryPlaceholder(output.valueKey)",
        "fieldSelectOptions(field)",
        ":placeholder=\"dataEntryPlaceholder(field.key)\"",
        ":placeholder=\"dataEntryPlaceholder(output.valueKey)\"",
        "Object.assign(this.dataEntryValues, this.blankDataEntryValues(), parsedValues)",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_data_entry_table_updates_patch_file_text_instead_of_regenerating_it():
    source = _source()
    select_tab_method = re.search(
        r"selectDataEntryTab\(tabId\) \{(?P<body>.*?)\n    \},",
        source,
        re.DOTALL,
    )
    table_input_method = re.search(
        r"onTableFieldInput\(fieldKey, value\) \{(?P<body>.*?)\n    \},",
        source,
        re.DOTALL,
    )

    assert select_tab_method, "selectDataEntryTab method is missing"
    assert table_input_method, "onTableFieldInput method is missing"
    assert "syncGeneratedInpcContent" not in select_tab_method.group("body")
    assert "patchInpcContentFromTableValues(fieldKey ? [fieldKey] : [])" in table_input_method.group("body")
    assert 'from "@/utils/epolyscat-input-script";' in source
    assert "dataEntryValues: {\n      deep: true" not in source


def test_uploaded_input_file_contents_seed_data_entry_table():
    source = _source()
    input_script_source = (
        Path(__file__).resolve().parents[1]
        / "src"
        / "utils"
        / "epolyscat-input-script.js"
    ).read_text()

    expected_hooks = [
        'import { InputService } from "@/service/epolyscat-service";',
        "parseEPolyScatInputScript as parseEPolyScatInputScriptFromContents",
        "async onInputStorageFilesSelected(files)",
        "async onInputComputerFilesSelected(event)",
        "async addFilesToActiveInput(files)",
        "readInputFileContents(file)",
        "loadInpcContentFromActiveInputFile",
        "loadInpcContentFromFile(file)",
        "InputService.fetchFileContents(readableFile)",
        "data-product-uri",
        "parseEPolyScatInputScriptFromContents(",
    ]

    for hook in expected_hooks:
        assert hook in source

    for hook in [
            "parseFileNameRecord",
            "parseConvertRecord",
            "parseAsyPolNode",
            "parseEngFormNode",
            "parseEPolyScatDocument",
            "serializeEPolyScatDocument",
            "normalizeEPolyScatRecord",
        "tokenizeEPolyScatLine",
        "stripEPolyScatComment",
    ]:
        assert hook in input_script_source


def test_uploaded_computer_files_keep_blob_identity_for_file_reader():
    source = _source()
    add_files_method = re.search(
        r"async addFilesToActiveInput\(files\) \{(?P<body>.*?)\n    \},\n    removeInputFile",
        source,
        re.DOTALL,
    )
    service_source = (
        Path(__file__).resolve().parents[1]
        / "src"
        / "service"
        / "epolyscat-service.js"
    ).read_text()

    assert add_files_method, "addFilesToActiveInput method is missing"
    assert "file.isFromComputer ? file : {" in add_files_method.group("body")
    assert "Object.assign(preparedFile, {" in add_files_method.group("body")
    assert "file instanceof Blob" in service_source
    assert "file.downloadURL" in service_source


def test_uploaded_computer_file_reading_preserves_blob_identity():
    source = _source()
    read_method = re.search(
        r"async readInputFileContents\(file\) \{(?P<body>.*?)\n    \},\n    async loadInpcContentFromActiveInputFile",
        source,
        re.DOTALL,
    )

    assert read_method, "readInputFileContents method is missing"
    body = read_method.group("body")
    assert "file.isFromComputer ? file : {" in body
    assert "...file" in body
    assert "InputService.fetchFileContents(readableFile)" in body


def test_storage_file_reading_accepts_user_storage_data_product_uri_key():
    source = _source()
    read_method = re.search(
        r"async readInputFileContents\(file\) \{(?P<body>.*?)\n    \},\n    async loadInpcContentFromActiveInputFile",
        source,
        re.DOTALL,
    )
    service_source = (
        Path(__file__).resolve().parents[1]
        / "src"
        / "service"
        / "epolyscat-service.js"
    ).read_text()

    assert read_method, "readInputFileContents method is missing"
    assert 'file["data-product-uri"]' in read_method.group("body")
    assert 'file["data-product-uri"]' in service_source
    assert "const dataProductURI" in service_source
    assert "const downloadURL" in service_source
    assert "encodeURIComponent(dataProductURI)" in service_source
    assert "if (!downloadURL)" in service_source


def test_storage_inp_files_download_even_when_mime_is_octet_stream():
    service_source = (
        Path(__file__).resolve().parents[1]
        / "src"
        / "service"
        / "epolyscat-service.js"
    ).read_text()

    assert "isKnownPlaintext" in service_source
    assert 'fileMetadata.isPlaintext(file.name) === true' in service_source
    assert "plaintextFileExtensions" in service_source
    assert "/\\.(inp|in|log|out|stdout|stderr|slurm)$/i" in service_source
    assert 'mimeType.startsWith("text")' in service_source
    assert "isKnownPlaintext || !mimeType" in service_source


def test_compute_resource_options_are_filtered_by_allocation_deployments_and_named():
    source = _resource_settings_source()

    expected_hooks = [
        "ApplicationDeploymentService.list",
        "groupResourceProfileId: this.groupResourceProfileId",
        "deploymentIds",
        "includedComputeResourceIds = new Set(deploymentIds)",
        "getComputeResourceDisplayName",
        "computeResourceName.hostName",
        "computeResourceName.name",
        "computeResourceName.value",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "if (deploymentIds.length > 0)" not in source


def test_create_run_reuses_run_store_save_and_submit_flow():
    source = _source()

    expected_hooks = [
        'this.$store.dispatch("run/createRun"',
        'this.$store.dispatch("run/submitRun"',
        "viewIds: this.viewIds",
        "router.push(`/runs/${this.run.parentRunId || this.workflowParentRunId || this.run.id}`)",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "RunService.createRun" not in source
    assert "directedit" not in source
    assert "inputTable" not in source


def test_resource_controls_are_searchable_input_dropdowns():
    source = _resource_settings_source()

    expected_hooks = [
        'class="resource-search-dropdown"',
        'class="resource-search-input"',
        'v-on:focus="openResourceMenu(\'allocation\')"',
        'v-on:focus="openResourceMenu(\'compute\')"',
        'v-on:input="onAllocationSearchInput"',
        'v-on:input="onComputeResourceSearchInput"',
        "showAllocationMenu",
        "showComputeResourceMenu",
        "allocationSearch",
        "computeResourceSearch",
        "filteredAllocationOptions",
        "filteredComputeResourceOptions",
        "selectAllocationOption",
        "selectComputeResourceOption",
        "loadResourceOptions",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert 'id="allocation"\n                  v-model="groupResourceProfileId"' not in source
    assert 'id="compute-resource"\n                  v-model="computeResourceId"' not in source
    assert "adpf-group-resource-profile-selector" not in source
    assert "adpf-experiment-compute-resource-selector" not in source


def test_resource_dropdowns_require_selecting_available_options():
    source = _resource_settings_source()

    forbidden_hooks = [
        "canUseTypedAllocation",
        "canUseTypedComputeResource",
        "applyTypedAllocationSearch",
        "applyTypedComputeResourceSearch",
        "typed-resource-option",
        "No matches",
    ]

    for hook in forbidden_hooks:
        assert hook not in source

    assert "Select an available option" in source


def test_resource_dropdown_options_do_not_render_identifier_meta():
    source = _resource_settings_source()

    assert "resource-option-meta" not in source
    assert '{{ option.value }}</span>' not in source


def test_resource_dropdown_valid_icon_is_offset_from_toggle_button():
    source = _resource_settings_source()

    assert ".resource-search-input.is-valid" in source
    assert "background-position: right 38px center;" in source
    assert "padding-right: 62px;" in source


def test_resource_area_uses_figma_width_instead_of_narrow_scaffold():
    source = _resource_settings_source()

    expected_hooks = [
        ".resource-settings-grid",
        "max-width: 1017px;",
        "grid-template-columns: minmax(0, 438px) minmax(0, 438px);",
        "width: calc(100% - 86px);",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "max-width: 860px;" not in source
    assert "width: 720px;" not in source


def test_resource_queue_card_uses_the_wider_card_internally():
    source = _resource_settings_source()

    expected_hooks = [
        "queue-card-header",
        "grid-template-columns: minmax(0, 1fr) auto;",
        ".queue-selector",
        "position: relative;",
        "grid-template-columns: minmax(92px, 1fr) 32px minmax(92px, 1fr) minmax(118px, 1fr) minmax(176px, 1.35fr);",
        "column-gap: 28px;",
        ".queue-input",
        "width: 100%;",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "grid-template-columns: 70px 32px 74px 88px 142px;" not in source
    assert "position: absolute;\n  right: 16px;\n  top: 8px;" not in source


def test_data_entry_output_and_script_long_values_can_wrap_or_scroll():
    source = _source()

    expected_hooks = [
        ".output-record-description",
        ".data-entry-script-preview",
        "overflow-wrap: anywhere;",
        "overflow: auto;",
        "white-space: pre;",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_compute_resources_are_empty_until_allocation_is_selected():
    source = _resource_settings_source()

    assert "if (!this.groupResourceProfileId)" in source
    assert "this.computeResourceOptions = [];" in source
    assert "Select allocation first" in source


def test_queue_settings_are_loaded_from_selected_application_deployment():
    source = _resource_settings_source()

    expected_hooks = [
        "ApplicationDeploymentService.getQueues",
        "lookup: applicationDeployment.appDeploymentId",
        "applicationDeployments.find",
        "deployment.computeHostId === this.computeResourceId",
        "queueOptions",
        "defaultQueue",
        "setDefaultQueue",
        "normalizeQueueOptions",
        "isQueueAllowedForSelection",
        "getBatchQueueResourcePolicy",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "Queue-skx-normal" not in source


def test_compute_resource_search_clears_stale_queue_options():
    source = _resource_settings_source()
    method = re.search(
        r"onComputeResourceSearchInput\(\) \{(?P<body>.*?)\n    \},",
        source,
        re.DOTALL,
    )
    assert method, "onComputeResourceSearchInput method is missing"
    body = method.group("body")

    assert "this.computeResourceId = null;" in body
    assert "this.clearQueueOptions();" in body


def test_queue_options_are_preserved_when_selected_resource_has_no_default_queue():
    source = _resource_settings_source()
    method = re.search(
        r"setDefaultQueue\(\) \{(?P<body>.*?)\n    \},",
        source,
        re.DOTALL,
    )
    assert method, "setDefaultQueue method is missing"
    body = method.group("body")

    assert "clearQueueSelection" in source
    assert "this.clearQueueSelection();" in body
    assert "this.clearQueueOptions();" not in body


def test_queue_switch_button_and_core_node_lock_follow_old_queue_editor():
    source = _resource_settings_source()

    expected_hooks = [
        'v-on:click="toggleQueueMenu"',
        'class="queue-menu"',
        "showQueueMenu",
        "filteredQueueOptions",
        "selectQueueOption",
        "queueResourceLockEnabled",
        "toggleQueueResourceLock",
        "onNodeCountInput",
        "onCoreCountInput",
        "selectedQueueDefault.cpuPerNode",
        "nodeCount * this.selectedQueueDefault.cpuPerNode",
        "Math.ceil(coreCount / this.selectedQueueDefault.cpuPerNode)",
        "b-icon :icon=\"queueResourceLockEnabled ? 'lock' : 'unlock'\"",
        ".queue-lock-button .b-icon",
        "height: 16px;",
        "width: 16px;",
    ]

    for hook in expected_hooks:
        assert hook in source

    assert "arrow-repeat" not in source
    assert "selectNextQueue" not in source


def test_create_run_preserves_required_validation_messages():
    source = _source() + _resource_settings_source()

    expected_messages = [
        "Run root is required",
        "Run input is required",
        "Compute resource is required",
        "Group resource profile is required",
        "A queue must be selected",
        "Core count must be greater than zero",
        "Node count must greater than zero",
        "Wall time limit must be specified",
        "Total physical memory must be specified",
    ]

    for message in expected_messages:
        assert message in source


def test_create_run_view_breadcrumb_links_back_to_view_runs():
    source = _source()

    assert '<b-breadcrumb-item to="/views">Views</b-breadcrumb-item>' in source
    assert ':to="`/runs?viewId=${view.viewId}`"' in source


def test_clone_run_preserves_group_resource_profile():
    source = Path(
        CREATE_RUN.parents[2] / "store" / "modules" / "run-storage.store.js"
    ).read_text()

    method = re.search(
        r"getStrippedRun: \([^)]*\) => \{(?P<body>.*?)\n    \}",
        source,
        re.DOTALL,
    )
    assert method, "getStrippedRun getter is missing"
    body = method.group("body")

    assert "groupResourceProfileId: run.groupResourceProfileId" in body
    assert body.index("groupResourceProfileId") < body.index("computeResourceId")


def test_routes_have_unique_new_run_related_names():
    runs_new_route = _route_block("/runs/new")
    create_run_alias_route = _route_block("/create-run")
    view_runs_route = _route_block("/views/:viewId")
    tutorial_runs_route = _route_block("/tutorials")

    assert "name: 'CreateRun'" in runs_new_route
    assert "name: 'CreateRunAlias'" in create_run_alias_route
    assert "name: 'ViewRuns'" in view_runs_route
    assert "name: 'TutorialRuns'" in tutorial_runs_route


def test_runs_new_route_loads_create_run_page():
    route = _route_block("/runs/new")

    assert "name: 'CreateRun'" in route
    assert "component: CreateRun" in route


def test_legacy_new_run_route_loads_old_run_page_before_dynamic_run_route():
    source = _router_source()
    route = _route_block("/runs/new/legacy")

    assert "name: 'LegacyCreateRun'" in route
    assert "component: Run" in route
    assert source.index("path: '/runs/new/legacy'") < source.index("path: '/runs/:runId'")


def test_run_detail_route_still_loads_run_page():
    route = _route_block("/runs/:runId")

    assert "name: 'Run'" in route
    assert "component: ViewRun" in route
