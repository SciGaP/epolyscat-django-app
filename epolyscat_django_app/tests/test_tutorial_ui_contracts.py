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


RUNS = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "components"
    / "Pages"
    / "Runs.vue"
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


WORKFLOW_FILE_LINKING = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "utils"
    / "workflow-file-linking.js"
)


TUTORIAL_CLONE = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "utils"
    / "tutorial-clone.js"
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


APP_LEFT_NAV = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "components"
    / "blocks"
    / "AppLeftNav.vue"
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


def _workflow_file_linking_source():
    return WORKFLOW_FILE_LINKING.read_text()


def _tutorial_clone_source():
    return TUTORIAL_CLONE.read_text()


def _input_store_source():
    return INPUT_STORE.read_text()


def _frontend_store_source():
    return FRONTEND_STORE.read_text()


def _run_store_source():
    return RUN_STORE.read_text()


def _view_store_source():
    return VIEW_STORE.read_text()


def _app_left_nav_source():
    return APP_LEFT_NAV.read_text()


def _view_run_source():
    return VIEW_RUN.read_text()


def _runs_source():
    return RUNS.read_text()


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


def test_tutorial_route_loads_real_shared_view_and_hides_mutations():
    service_source = _epolyscat_service_source()
    fetch_view = re.search(
        r"async fetchView\(\{viewId\}\) \{(?P<body>.*?)\n    \},",
        service_source,
        re.DOTALL,
    )
    assert fetch_view, "ViewService.fetchView is missing"
    fetch_view_body = fetch_view.group("body")

    assert "views/tutorials/" in fetch_view_body
    assert "runs/tutorial_runs" not in fetch_view_body
    assert "id: -1" not in fetch_view_body
    assert 'v-if="view != null && !isTutorials"' in _runs_source()


def test_tutorial_navigation_and_view_store_keep_shared_view_isolated():
    nav_source = _app_left_nav_source()
    view_store_source = _view_store_source()

    assert 'router-link :to="`/tutorials`"' in nav_source
    assert "state.viewListMap[queryString]" in view_store_source
    assert "viewIds.map" in view_store_source


def test_tutorial_list_uses_tutorial_format_status_and_clone_columns():
    runs_source = _runs_source()
    service_source = _epolyscat_service_source()

    assert ':columns="runColumns"' in runs_source
    assert "runColumns()" in runs_source
    assert "runSorters()" in runs_source
    assert "['name', 'Tutorial']" in runs_source
    assert "['inputFormat', 'Input Format']" in runs_source
    assert "['status', 'Status']" in runs_source
    assert "['actions', 'Actions']" in runs_source
    assert 'v-if="isTutorials && item.description"' in runs_source
    assert 'class="tutorial-name-cell"' in runs_source
    assert ".tutorial-name-cell" in runs_source
    assert "flex-direction: column" in runs_source
    assert "item.inputFormat" in runs_source
    assert "tutorialInputFormat" in service_source


def test_tutorial_list_failure_never_falls_back_to_user_runs_and_can_retry():
    runs_source = _runs_source()

    assert "tutorialLoadError" in runs_source
    assert "retryTutorialLoad" in runs_source
    assert "if (this.isTutorials && this.view == null) return [];" in runs_source
    assert 'v-if="isTutorials && tutorialLoadError"' in runs_source
    assert "Tutorials could not be loaded." in runs_source


def test_tutorial_detail_preserves_manifest_file_groups():
    source = _view_run_source()

    assert "if (this.run && this.run.isTutorial)" in source
    assert "return this.presentation.file_groups" in source
    assert "tutorial-file-groups" in source
    assert "<h2>{{ group.label }}</h2>" in source
    assert 'v-for="file in group.files"' in source
    assert 'label: "Reference"' not in source


def test_tutorial_detail_offers_clone_without_edit_or_delete_controls():
    source = _view_run_source()

    assert "Clone Tutorial" in source
    assert 'v-if="run && run.isTutorial"' in source
    assert "cloneTutorial" in source
    assert "`/runs/new?clonedFrom=${this.run.id}`" in source
    assert (
        "run && run.isTutorial ? &quot;/tutorials&quot; : &quot;/runs&quot;"
        in source
    )
    assert 'run && run.isTutorial ? "Tutorials" : "Experiments"' in source
    assert source.count("<RunResource") == 1


def test_tutorial_clone_hydrates_packaged_inputs_for_normal_run_creation():
    create_source = _source()
    run_store_source = _run_store_source()
    linking_source = _workflow_file_linking_source()
    tutorial_clone_source = _tutorial_clone_source()

    assert "cloningTutorial" in create_source
    assert "applyCloneSourceRun" in create_source
    assert "hydrateTutorialCloneInputs" in create_source
    assert "prepareTutorialCloneInputs" in create_source
    assert "RunService.fetchViewableContent" in create_source
    assert "dataProductURI: null" in tutorial_clone_source
    assert "downloadURL: null" in tutorial_clone_source
    assert "unchanged: false" in tutorial_clone_source
    assert "tutorialStandalone" in create_source
    assert "tutorialRuntimeDataName" in tutorial_clone_source
    assert "tutorialConvertFormat" in tutorial_clone_source
    assert "tutorialCloneHydrating" in create_source
    assert "tutorialCloneHydrationError" in create_source
    assert "retryTutorialCloneHydration" in create_source
    assert "tutorialCloneBlocked" in create_source
    assert "patchEPolyScatInputScript" in tutorial_clone_source
    assert 'changedKeys: ["convertSource", "convertFormat"]' in tutorial_clone_source
    assert "this.applyInpcContentToDataEntryValues(hydratedClone.scriptContent)" in create_source
    assert 'file.name !== "ePolyScat_Input_Data"' not in create_source
    assert "convertFormat = \"\"" in linking_source
    assert "convertFormat ||" in linking_source

    for field in (
        "runMode",
        "moduleApplication",
        "workflowStage",
        "workflowApplication",
        "utilityApplication",
        "workflowMetadata",
    ):
        assert f"{field}: run.{field}" in run_store_source


def test_integrated_clone_keeps_application_identity_and_both_save_guards():
    source = CREATE_RUN.read_text()
    helper = re.search(r'applyCloneSourceRun\(sourceRun\) \{(?P<body>.*?)\n    \},', source, re.DOTALL).group('body')
    assert 'applicationSelectionForRun(sourceRun)' in helper
    assert 'sourceRun.utilityApplication || "CnvMath"' not in helper
    assert 'processing || !applicationSelectionAvailable || tutorialCloneBlocked' in source
    assert 'if (!this.applicationSelectionAvailable || this.tutorialCloneBlocked) return;' in source
    runs = RUNS.read_text()
    assert '<ListView class="portal-list-table"' in runs
    assert 'run.experimentId == this.experimentId' in runs
    assert 'this.isTutorials ? "Tutorials"' in runs
