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
        r"<b-button[^>]*@click=\"goToNewRun\"[^>]*>\s*New Run\s*</b-button>",
        source,
        re.DOTALL,
    )

    assert "create-run-type-modal" not in source
    assert "chooseWorkflowRun" not in source
