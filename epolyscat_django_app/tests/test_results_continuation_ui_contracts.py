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


INPUT_STORE = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "store"
    / "modules"
    / "input-storage.store.js"
)


INPUT_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "utils"
    / "epolyscat-input-script.js"
)


FRONTEND_STORE = Path(__file__).resolve().parents[1] / "src" / "store" / "index.js"


RUN_STORE = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "store"
    / "modules"
    / "run-storage.store.js"
)


VIEW_STORE = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "store"
    / "modules"
    / "view-storage.store.js"
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


def _input_store_source():
    return INPUT_STORE.read_text()


def _frontend_store_source():
    return FRONTEND_STORE.read_text()


def _run_store_source():
    return RUN_STORE.read_text()


def _view_store_source():
    return VIEW_STORE.read_text()


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



def test_create_run_and_workflow_pages_offer_run_type_switches():
    create_source = _source()
    workflow_source = _workflow_source()

    assert "run-type-switcher-row" in create_source
    assert "run-selection-column run-type-column" in create_source
    assert "selectRunType(runType.id)" in create_source
    assert "Switch to Workflow Run" not in create_source
    assert "switchToModuleRun" in workflow_source
    assert "Switch to Module Run" in workflow_source


def test_dedicated_application_ids_do_not_fall_back_to_epolyscat_deployments():
    source = _source()
    resource_source = _resource_settings_source()
    workflow_source = _workflow_source()
    computed = re.search(
        r"resourceApplicationModuleId\(\) \{(?P<body>.*?)\n    \},",
        source,
        re.DOTALL,
    )
    assert computed, "resourceApplicationModuleId computed property is missing"
    body = computed.group("body")

    assert 'return this.$store.getters["settings/gaussian16ApplicationModuleId"];' in body
    assert 'return this.$store.getters["settings/openmolcasApplicationModuleId"];' in body
    assert "|| epolyscatApplicationModuleId" not in body
    assert 'return this.applicationModuleId;' in resource_source
    assert 'settings/epolyscatApplicationModuleId' not in resource_source
    assert ':application-module-id="epolyscatApplicationModuleId"' in workflow_source


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
        "this.run.presentation",
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

    assert "PlotService.getViewables" not in source


def test_view_run_parent_workflow_uses_active_child_for_file_catalog():
    source = _view_run_source()

    expected_hooks = [
        "fileCatalogRunId",
        'this.presentation.mode === "workflow"',
        "!this.isWorkflowChildRun",
        "this.activeWorkflowStage.child_run_id",
        "PlotService.getInputFiles({ runId: this.fileCatalogRunId })",
        "InputService.fetchOutputs(this.fileCatalogRunId)",
        "runId: this.fileCatalogRunId",
        "runIds: [this.fileCatalogRunId]",
    ]

    for hook in expected_hooks:
        assert hook in source


def test_view_run_reuses_run_presentation_and_output_catalog():
    source = _view_run_source()

    assert "this.run.presentation" in source
    assert "RunService.fetchRunPresentation" not in source
    assert "PlotService.getViewables" not in source
    assert source.index("await this.loadRunFileCatalog()") < source.index(
        "await this.loadWorkflowContinuation()"
    )


def test_frontend_stores_coalesce_concurrent_list_requests():
    run_source = _run_store_source()
    view_source = _view_store_source()

    assert "let pendingRunsRequest = null" in run_source
    assert "if (pendingRunsRequest)" in run_source
    assert "pendingRunsRequest = null" in run_source

    assert "const pendingViewRequests = {}" in view_source
    assert "pendingViewRequests[queryString]" in view_source
    assert "delete pendingViewRequests[queryString]" in view_source


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
        "workflow-continuation-preview",
        "workflowContinuation.nextStagePreview",
        "continuationPreviewOutputName",
        'label: "Continue to"',
        'label: "Inherited input"',
    ]
    expected_service_hooks = [
        "fetchWorkflowContinuation",
        "continueWorkflow",
        "workflow_continuation/",
        "workflow_parent_run_id",
        "next_child_run_id",
        "source_run_id",
        "next_stage_preview",
        "nextStagePreview",
    ]

    for hook in expected_view_hooks:
        assert hook in source
    for hook in expected_service_hooks:
        assert hook in service_source


def test_scientific_verification_is_visible_and_blocks_silent_output_inheritance():
    create_run_source = _source()
    view_run_source = _view_run_source()
    service_source = _epolyscat_service_source()

    for hook in [
        "workflowOutputBindingMessage",
        'variant="warning"',
        'backendBinding.status === "blocked"',
        "backendBinding.message",
    ]:
        assert hook in create_run_source

    for hook in [
        'label: "Science"',
        "workflowContinuation.scientificVerification",
        "verification.message",
    ]:
        assert hook in view_run_source

    for hook in [
        "scientificVerification: data.scientific_verification",
        "provenance: data.provenance",
    ]:
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
        "file.plottable === true",
        "file.plot_contract",
        "applyPlotContractForFile",
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
