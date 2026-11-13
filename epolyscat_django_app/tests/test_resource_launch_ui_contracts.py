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


def test_hidden_queue_editor_cannot_overwrite_visible_resource_values():
    source = _resource_settings_source()

    assert '<adpf-queue-settings-editor' in source
    assert 'v-on:input="onQueueSettingEditor"' not in source
    assert "onQueueSettingEditor(evt)" not in source
    assert "applyQueueDefaults(queueDefault)" in source
