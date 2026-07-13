"""Domain rules for continuing completed runs through the ePolyScat workflow."""


WORKFLOW_STAGES = ("Data_Gen", "ePolyScat_Run", "Analysis")

MODULE_STAGES = {
    "gaussian16": ("Gaussian16", "Data_Gen"),
    "openmolcas": ("OpenMolcas", "Data_Gen"),
    "epolyscat": ("ePolyScat", "ePolyScat_Run"),
}

UTILITY_APPLICATIONS = {
    "cnvmath": "CnvMath",
    "cnvmatlab": "CnvMatLab",
    "cnvlinfull": "CnvLinFull",
    "moldenmerge": "MoldenMerge",
    "nrfpad": "NRFPAD",
    "cube2igor": "Cube2igor",
}

STAGE_ALIASES = {
    "data-generation": "Data_Gen",
    "data_generation": "Data_Gen",
    "data_gen": "Data_Gen",
    "bound": "ePolyScat_Run",
    "epolyscat-dmat": "ePolyScat_Run",
    "epolyscat_run": "ePolyScat_Run",
    "analysis": "Analysis",
}


def _selector_inputs(run):
    if not getattr(run, "pk", None):
        return {}
    try:
        return {
            input_instance.name: input_instance.value
            for input_instance in run.inputs.all()
        }
    except (AttributeError, ValueError):
        return {}


def _canonical_stage(value):
    value = str(value or "")
    return STAGE_ALIASES.get(value.lower(), value if value in WORKFLOW_STAGES else "")


def _next_stage(stage):
    try:
        return WORKFLOW_STAGES[WORKFLOW_STAGES.index(stage) + 1]
    except (ValueError, IndexError):
        return None


def classify_run(run):
    inputs = _selector_inputs(run)
    run_mode = str(getattr(run, "run_mode", "") or "").lower()
    has_structured_selector = any(
        getattr(run, field, "")
        for field in (
            "module_application",
            "workflow_stage",
            "workflow_application",
            "utility_application",
        )
    )
    if not has_structured_selector and inputs.get("Calculation_Type"):
        run_mode = str(inputs["Calculation_Type"]).lower()

    if run_mode == "workflow":
        metadata = getattr(run, "workflow_metadata", None) or {}
        stage = _canonical_stage(
            getattr(run, "workflow_stage", "")
            or inputs.get("Application_Workflow")
        )
        if metadata.get("isWorkflowPlan") is True or not stage:
            return None
        if stage == "Data_Gen":
            application = (
                getattr(run, "workflow_application", "")
                or inputs.get("Data_Gen")
                or inputs.get("Workflow_Application")
            )
        elif stage == "Analysis":
            application = (
                getattr(run, "utility_application", "")
                or inputs.get("Application_Utility")
            )
        else:
            application = getattr(run, "module_application", "") or "ePolyScat"
        if not application:
            return None
        return {
            "source_stage": stage,
            "source_application": application,
            "next_stage": _next_stage(stage),
        }

    if run_mode == "utility":
        application = (
            getattr(run, "utility_application", "")
            or inputs.get("Application_Utility")
        )
        canonical_application = UTILITY_APPLICATIONS.get(str(application).lower())
        if not canonical_application:
            return None
        return {
            "source_stage": "Analysis",
            "source_application": canonical_application,
            "next_stage": None,
        }

    application = (
        getattr(run, "module_application", "")
        or inputs.get("EPOLYSCAT_Application_Module")
    )
    module_stage = MODULE_STAGES.get(str(application).lower())
    if not module_stage:
        return None
    canonical_application, stage = module_stage
    return {
        "source_stage": stage,
        "source_application": canonical_application,
        "next_stage": _next_stage(stage),
    }


def continuation_eligibility(run, request):
    classification = classify_run(run)

    def unavailable(reason, message):
        return {
            "eligible": False,
            "reason": reason,
            "message": message,
            **(classification or {}),
        }

    if getattr(run, "deleted", False):
        return unavailable("deleted_run", "Deleted runs cannot start a workflow.")
    if classification is None:
        return unavailable(
            "unsupported_run_type",
            "This run type is not part of the ePolyScat workflow.",
        )
    if classification["next_stage"] is None:
        return unavailable(
            "terminal_stage",
            "This run is already at the final implemented workflow stage.",
        )

    latest_execution = getattr(run, "latest_execution", None)
    if latest_execution is None:
        return unavailable(
            "no_execution",
            "Submit and complete this run before continuing it in a workflow.",
        )

    try:
        status = latest_execution.get_airavata_experiment_status(request)
    except Exception:
        return unavailable(
            "status_unavailable",
            "The source run status could not be confirmed.",
        )
    if str(status).upper() != "COMPLETED":
        return unavailable(
            "run_not_completed",
            "Only completed runs can continue into the next workflow stage.",
        )

    return {
        "eligible": True,
        "reason": "",
        "message": "",
        **classification,
    }
