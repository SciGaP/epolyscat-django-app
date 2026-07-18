import json
import logging
import os
from io import StringIO
from urllib.parse import urlencode

from .airavata_grpc import Project, django_user, experiment_util, user_storage
from django.apps import apps
from django.db import transaction
from django.db.models import Q
from django.conf import settings
from django.utils.text import get_valid_filename
from rest_framework import reverse, serializers, validators

from epolyscat_django_app import models

logger = logging.getLogger(__name__)

class FileSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.File
        fields = ['name', 'data_product_uri']

class InputSerializer(serializers.ModelSerializer):
    files = serializers.SerializerMethodField()

    class Meta:
        model = models.Input
        fields = ['type', 'name', 'value', 'files']

    def get_files(self, input_instance):
        files = models.File.objects.filter(input=input_instance)

        return FileSerializer(files, many=True).data


class UniqueToUserValidator(validators.UniqueValidator):
    requires_context = True

    def __init__(self, queryset, user_field, message=None, lookup="exact"):
        self.user_field = user_field
        super().__init__(queryset, message=message, lookup=lookup)

    def __call__(self, value, serializer_field):
        self.user = django_user(serializer_field.context["request"])
        return super().__call__(value, serializer_field)

    def filter_queryset(self, value, queryset, field_name):
        # filter by current user
        queryset = queryset.filter(**{self.user_field: self.user})
        return super().filter_queryset(value, queryset, field_name)


class ExperimentIdRelatedField(serializers.PrimaryKeyRelatedField):
    def get_queryset(self):
        request = self.context["request"]
        return models.Experiment.objects.filter(owner=django_user(request), deleted=False)


class RunSerializer(serializers.ModelSerializer):
    inputs = serializers.SerializerMethodField()
    executions = serializers.SlugRelatedField(
        slug_field="airavata_experiment_id", read_only=True, many=True
    )
    experiment = ExperimentIdRelatedField(required=False, allow_null=True)

    #root = serializers.CharField(max_length=100, required=False)
    #directedit = serializers.CharField(
    #    style={"base_template": "textarea.html"},
    #    allow_blank=True,
    #    write_only=True,
    #    required=False,
    #)
    #inpc_download_url = serializers.SerializerMethodField()
    status = serializers.SerializerMethodField()
    job_status = serializers.SerializerMethodField()
    is_tutorial = serializers.SerializerMethodField()
    job_id = serializers.SerializerMethodField()
    resource = serializers.SerializerMethodField()
    presentation = serializers.SerializerMethodField()
    #resource_short = serializers.SerializerMethodField()
    #executions = serializers.SlugRelatedField(
    #    slug_field="airavata_experiment_id", read_only=True, many=True
    #)
    #input_table = serializers.JSONField(allow_null=True, required=False)
    #can_resubmit = serializers.SerializerMethodField()
    #cancelable = serializers.SerializerMethodField()

    class Meta:
        model = models.Run
        fields = (
            "id", "owner", "name", "description", "airavata_project_id", "views",
            "created", "updated", "deleted", 'is_email_notification_on',
            "group_resource_profile_id", "compute_resource_id",
            "queue_name", "core_count", "node_count", "walltime_limit", "total_physical_memory",
            "run_mode", "module_application", "workflow_stage", "workflow_application",
            "utility_application", "workflow_metadata",
            "parent_run", "workflow_source_run", "workflow_step_order", "workflow_step_status",
            "inputs", "executions", "status", "job_status", "is_tutorial", "job_id",
            "resource", "presentation", "experiment",
            #"directedit", "inpc_download_url", "cancelable","can_resubmit", "input_table", "root",
            #"number", "root", "resource", #"resource_short", "job_id",
        )
        #read_only_fields = ("deleted", "number", "experiment", "name")

    #def to_representation(self, instance):
    #    rep = super().to_representation(instance)
    #    rep["root"] = instance.root.root
    #    if instance.input_table is not None:
    #        rep["input_table"] = json.loads(instance.input_table)
    #    return rep

    @transaction.atomic
    def create(self, validated_data):
        request = self.context["request"]

        projects = request.airavata_client.getUserProjects(
            request.authz_token,
            settings.GATEWAY_ID,
            request.user.username,
            -1,
            0
        )
        epolyscat_project_choices = ([p for p in projects if "EPOLYSCAT_app_project" in p.projectID] or
            [p for p in projects if "Default_Project" in p.projectID] or
            [p for p in projects if "Default" in p.projectID] or
            [p for p in projects if "default" in p.projectID])

        if len(epolyscat_project_choices) > 0:
            airavata_project_id = epolyscat_project_choices[0].projectID
        else:
            new_project = Project(
                owner=request.user.username,
                gatewayId=settings.GATEWAY_ID,
                name="EPOLYSCAT app project",
            )

            airavata_project_id = request.airavata_client.createProject(
                request.authz_token,
                settings.GATEWAY_ID,
                new_project
            )

        # print(projects, [p for p in projects if p.projectID.startswith("Default_Project")])
        # default_project = next(p for p in projects if p.projectID.startswith("Default_Project"))

        # airavata_project_id = default_project.projectID

        return models.Run.objects.create(
            **validated_data,
            airavata_project_id=airavata_project_id,
            directory=""
        )

        '''
        root = get_valid_filename(validated_data.pop("root"))
        runs_root, created = models.RunsRoot.objects.get_or_create(
            root=root, owner=request.user
        )
        if created:
            experiment = models.Experiment.objects.create(
                name=root, root=runs_root, owner=request.user
            )
            experiment.create_airavata_project(request)
            experiment.save()
        directedit = validated_data.pop("directedit", "")
        input_table = validated_data.pop("input_table", None)
        experiment = runs_root.experiment
        number = runs_root.get_next_run_number()
        run = models.Run.objects.create(
            **validated_data,
            root=runs_root,
            number=number,
            experiment=experiment,
            airavata_project_id=airavata_project_id,
            directory="",
        )
        run_dirs = ("Runs", runs_root.root, run.number)
        user_storage.create_user_dir(request, dir_names=run_dirs)
        # filepath is relative to user directory instead of the full path
        run.filepath = os.path.join(*run_dirs)
        if directedit.strip() != "":
            self._create_inpc_file(run, directedit)
        elif input_table is not None:
            self._create_inpc_file_input_table(run, input_table)
        run.save()
        return run

    @transaction.atomic
    def update(self, instance, validated_data):
        request = self.context["request"]
        view = self.context["view"]

        # Always update queue settings, even if resubmitting
        instance.queue_name = validated_data.get("queue_name", instance.queue_name)
        instance.core_count = validated_data.get("core_count", instance.core_count)
        instance.node_count = validated_data.get("node_count", instance.node_count)
        instance.walltime_limit = validated_data.get(
            "walltime_limit", instance.walltime_limit
        )
        instance.total_physical_memory = validated_data.get(
            "total_physical_memory", instance.total_physical_memory
        )

        if view.action not in ["resubmit"]:
            instance.group_resource_profile_id = validated_data.get(
                "group_resource_profile_id", instance.group_resource_profile_id
            )
            instance.compute_resource_id = validated_data.get(
                "compute_resource_id", instance.compute_resource_id
            )
            directedit = validated_data.pop("directedit", "")
            input_table = validated_data.pop("input_table", None)

            # if file exists, update it, else create it
            if instance.inpc_data_product_uri is not None and user_storage.exists(
                request, data_product_uri=instance.inpc_data_product_uri
            ):
                new_inpc_string = None

                # validation guarantees that one of 'directedit' or 'input_table' is available
                if directedit.strip() != "":
                    new_inpc_string = directedit
                elif input_table is not None:
                    new_inpc_string = self._create_inpc_string_from_input_table(
                        input_table
                    )
                    instance.input_table = json.dumps(input_table)

                if new_inpc_string is not None:
                    user_storage.update_data_product_content(
                        request,
                        data_product_uri=instance.inpc_data_product_uri,
                        fileContentText=new_inpc_string,
                    )
            else:
                if directedit.strip() != "":
                    self._create_inpc_file(instance, directedit)
                elif input_table is not None:
                    self._create_inpc_file_input_table(instance, input_table)
        return instance

    def validate(self, attrs):
        view = self.context["view"]
        if view.action == "submit":
            # Validate that execution parameters are provided
            # For now we won't worry about the correctness of the parameters,
            # just checking that they have a value
            submit_required_fields = [
                "group_resource_profile_id",
                "compute_resource_id",
                "queue_name",
                "core_count",
                "node_count",
                "walltime_limit",
                "total_physical_memory",
            ]
            for field in submit_required_fields:
                value = attrs.get(field, None)
                if value is None or value == "":
                    raise serializers.ValidationError(
                        f"{field} must be provided for submission"
                    )
        if view.action in ("create", "update"):
            directedit = attrs.get("directedit", "")
            input_table = attrs.get("input_table", None)
            directedit_provided = directedit is not None and directedit != ""
            input_table_provided = input_table is not None
            if directedit_provided and input_table_provided:
                raise serializers.ValidationError(
                    "Must not supply values for both directedit and input_table"
                )
            if not directedit_provided and not input_table_provided:
                raise serializers.ValidationError(
                    "Please provide one of 'directedit' or 'input_table' to specify the input file"
                )
        if view.action == "create":
            root = attrs.get("root", None)
            if root is None:
                raise serializers.ValidationError("'root' is required to create a run.")
        return attrs

    def get_inpc_download_url(self, instance):
        request = self.context["request"]
        if instance.inpc_data_product_uri is not None:
            return user_storage.get_download_url(
                request, data_product_uri=instance.inpc_data_product_uri
            )
        else:
            return None
    '''
    def get_inputs(self, run_instance: models.Run):
        inputs = models.Input.objects.filter(run=run_instance)

        return InputSerializer(inputs, many=True).data

    def get_presentation(self, run_instance: models.Run):
        mode = run_instance.run_mode or "module"
        metadata = run_instance.workflow_metadata or {}
        context_run = self._presentation_context_run(run_instance)
        schema = self._presentation_file_schema(run_instance, context_run)
        workflow_stage = schema["active_stage_id"]
        input_files = self._run_file_names(run_instance, "files")
        output_files = metadata.get("output_files") or schema["output_files"]
        target_states = {
            "input_columns": metadata.get("target_state_input_columns") or schema["input_columns"],
            "output_files": output_files,
        }
        selectable_inputs = input_files or schema["input_files"]
        selectable_files = selectable_inputs or output_files or [""]
        file_groups = [
            {"label": "Inputs", "files": selectable_inputs},
            {"label": "Outputs", "files": output_files + schema["output_extras"]},
        ]

        presentation = {
            "mode": mode,
            "subtitle": schema["subtitle"],
            "active_stage_id": workflow_stage,
            "target_states": target_states,
            "file_groups": file_groups,
            "selected_file": metadata.get("selected_file") or selectable_files[0],
            "code": metadata.get("code") or self._default_code(run_instance),
            "parameters": metadata.get("parameters") or [
                {"label": "Coupling mode", "value": "LS"},
                {"label": "Nuclear charge/Atomic number", "value": "2"},
                {"label": "Number of electrons", "value": "1"},
            ],
            "plot": metadata.get("plot") or {
                "file": schema["plot_file"],
                "x_axis": "0",
                "y_axis": "1",
                "flags": "-linY",
            },
            "plottable_file_names": self._plottable_file_names(),
            "applications": self._workflow_applications(run_instance.workflow_application),
            "stages": self._workflow_stages("epolyscat-dmat"),
        }

        if mode == "workflow":
            child_stages = self._workflow_child_stages(run_instance)
            presentation["stages"] = child_stages or self._workflow_stages(workflow_stage)
            presentation.update(self._workflow_progress(run_instance, presentation["stages"]))

        return presentation

    def _run_file_names(self, run_instance: models.Run, input_type):
        names = []
        inputs = models.Input.objects.filter(run=run_instance, type=input_type)
        for input_instance in inputs:
            names.extend(file.name for file in input_instance.files.all())
        return names

    def _plottable_file_names(self):
        # Kept in the presentation shape for older clients. New clients use the
        # per-output plottable and plot_contract fields from get_output_files.
        return []

    def _presentation_context_run(self, run_instance):
        if run_instance.run_mode != "workflow" or run_instance.parent_run_id is not None:
            return run_instance

        children = list(run_instance.workflow_steps.order_by("workflow_step_order", "id"))
        if not children:
            return run_instance

        metadata = run_instance.workflow_metadata or {}
        active_child_run_id = metadata.get("active_child_run_id")
        if active_child_run_id:
            active_child = next(
                (child for child in children if child.id == active_child_run_id),
                None,
            )
            if active_child is not None:
                return active_child

        active_child = next(
            (
                child for child in children
                if child.workflow_step_status in ("submitted", "running")
                or child.executions.exists()
            ),
            None,
        )
        return active_child or children[0]

    def _presentation_file_schema(self, display_run, context_run):
        application = self._presentation_application(context_run)
        stage = self._normalized_workflow_stage(context_run.workflow_stage)
        display_stage = context_run.workflow_stage or stage
        schema = self._file_schema_for_application(application)

        if display_run.run_mode == "workflow":
            subtitle = f"Workflow/{self._workflow_stage_label(display_stage)}"
            active_stage_id = display_stage
        elif display_run.run_mode == "utility":
            subtitle = f"Utilities/{application}"
            active_stage_id = application
        elif display_run.module_application:
            subtitle = f"Modules/{application}"
            active_stage_id = application
        else:
            subtitle = "Modules/EPOLYSCAT_DMAT"
            active_stage_id = "epolyscat-dmat"

        return {
            **schema,
            "subtitle": subtitle,
            "active_stage_id": active_stage_id,
        }

    def _presentation_application(self, run_instance):
        stage = self._normalized_workflow_stage(run_instance.workflow_stage)
        if run_instance.run_mode == "utility":
            return run_instance.utility_application or "CnvMath"
        if run_instance.run_mode == "workflow":
            if stage == "Data_Gen":
                return run_instance.workflow_application or "OpenMolcas"
            if stage == "Analysis":
                return run_instance.utility_application or "CnvMath"
            return run_instance.module_application or "ePolyScat"
        return run_instance.module_application or "ePolyScat"

    def _normalized_workflow_stage(self, stage):
        aliases = {
            "data-generation": "Data_Gen",
            "Data_Generation": "Data_Gen",
            "bound": "ePolyScat_Run",
            "epolyscat-dmat": "ePolyScat_Run",
            "analysis": "Analysis",
        }
        return aliases.get(stage, stage or "ePolyScat_Run")

    def _file_schema_for_application(self, application):
        epolyscat_inputs = [
            ["ns_001.c", "target", "ns_001.bsw"],
            ["nd_001.c", "knot.dat", "nd_001.bsw"],
            ["bound.nnn", "target.bsw", "mult_bnk", "pert_nnn.bsw"],
        ]
        schemas = {
            "Gaussian16": {
                "input_columns": [["Gaussian_Input"]],
                "input_files": ["Gaussian_Input"],
                "output_files": ["gaussian.log", "molden.dat"],
                "output_extras": [],
                "plot_file": "molden.dat",
            },
            "OpenMolcas": {
                "input_columns": [["Molcas_Input"]],
                "input_files": ["Molcas_Input"],
                "output_files": ["molcas.log", "molden.dat"],
                "output_extras": [],
                "plot_file": "molden.dat",
            },
            "ePolyScat": {
                "input_columns": epolyscat_inputs,
                "input_files": ["ePolyScat_Input_Data", "ePolyscat_Input_File"],
                "output_files": ["d.nnn", "zf_res", "ePolyScat_dmat.log"],
                "output_extras": ["Parameters", "bound_tab", "logs"],
                "plot_file": "cross_sections",
            },
            "MoldenMerge": {
                "input_columns": [["molden.dat"]],
                "input_files": ["molden.dat"],
                "output_files": ["merged_molden.dat"],
                "output_extras": [],
                "plot_file": "merged_molden.dat",
            },
        }
        utility_outputs = {
            "CnvMath": "mathematica_output.dat",
            "CnvMatLab": "matlab_output.dat",
            "CnvLinFull": "differential_cross_section.dat",
            "NRFPAD": "nrfpad.dat",
            "Cube2igor": "igor_plot.itx",
        }
        if application in utility_outputs:
            return {
                "input_columns": [],
                "input_files": [],
                "output_files": [utility_outputs[application]],
                "output_extras": [],
                "plot_file": utility_outputs[application],
            }
        return schemas.get(application, schemas["ePolyScat"])

    def _workflow_stage_label(self, stage):
        labels = {
            "Data_Gen": "Data Generation",
            "Data_Generation": "Data Generation",
            "data-generation": "Data Generation",
            "bound": "Bound",
            "ePolyScat_Run": "ePolyScat Run",
            "stgf": "STGF",
            "Analysis": "Post-processing",
            "analysis": "Post-processing",
            "Visualization": "Visualization",
            "visualization": "Visualization",
            "epolyscat-dmat": "EPOLYSCAT_DMAT",
        }
        return labels.get(stage, stage.replace("-", " ").title())

    def _workflow_applications(self, selected_application):
        selected = selected_application or "OpenMolcas"
        return [
            {
                "id": "Gaussian16",
                "label": "Gaussian16",
                "selected": selected == "Gaussian16",
                "parameters": [
                    {"label": "Input file", "name": "gaussian_input"},
                    {"label": "GPU_Version?", "name": "gpu_version", "value": "Yes"},
                ],
            },
            {
                "id": "OpenMolcas",
                "label": "OpenMolcas",
                "selected": selected == "OpenMolcas",
                "parameters": [
                    {"label": "Input file", "name": "openmolcas_input"},
                    {"label": "Optional files", "name": "openmolcas_optional"},
                    {"label": "Printing of orbitals", "name": "printing_orbitals", "value": "Yes"},
                ],
            },
        ]

    def _workflow_stages(self, active_stage_id):
        stage_aliases = {
            "data-generation": "Data_Gen",
            "Data_Generation": "Data_Gen",
            "bound": "ePolyScat_Run",
            "analysis": "Analysis",
            "visualization": "Visualization",
        }
        active_stage_id = stage_aliases.get(active_stage_id, active_stage_id)
        stages = [
            {"id": "Data_Gen", "label": "Data Generation"},
            {"id": "ePolyScat_Run", "label": "ePolyScat Run"},
            {"id": "Analysis", "label": "Post-processing"},
            {
                "id": "Visualization",
                "label": "Visualization",
                "local_only": True,
            },
        ]
        active_index = next(
            (index for index, stage in enumerate(stages) if stage["id"] == active_stage_id),
            1,
        )
        for index, stage in enumerate(stages):
            if index < active_index:
                stage["state"] = "complete"
            elif index == active_index:
                stage["state"] = "active"
            else:
                stage["state"] = "pending"
        return stages

    def _workflow_child_stages(self, run_instance):
        if run_instance.parent_run_id is not None:
            return self._workflow_child_stages(run_instance.parent_run)

        child_runs = list(run_instance.workflow_steps.order_by("workflow_step_order", "id"))
        if not child_runs:
            return []

        metadata = run_instance.workflow_metadata or {}
        imported_source = metadata.get("importedSource") or {}
        imported_stage = imported_source.get("stage")
        if run_instance.workflow_source_run_id and imported_stage:
            canonical_stages = ["Data_Gen", "ePolyScat_Run", "Analysis"]
            try:
                imported_index = canonical_stages.index(imported_stage)
            except ValueError:
                imported_index = -1
            if imported_index >= 0:
                stages = [
                    {
                        "id": stage,
                        "label": self._workflow_stage_label(stage),
                        "state": "not_included",
                        "status": "not_included",
                    }
                    for stage in canonical_stages[:imported_index]
                ]
                stages.append(
                    {
                        "id": imported_stage,
                        "label": self._workflow_stage_label(imported_stage),
                        "state": "complete",
                        "status": "imported",
                        "child_run_id": run_instance.workflow_source_run_id,
                        "application": imported_source.get("application") or "",
                        "imported": True,
                    }
                )
                stages.extend(self._workflow_children_to_stages(child_runs))
                return self._with_visualization_stage(stages)

        stages = self._workflow_children_to_stages(child_runs)
        return self._with_visualization_stage(stages)

    def _workflow_children_to_stages(self, child_runs):
        stages = []
        analysis_children = [
            child for child in child_runs if child.workflow_stage == "Analysis"
        ]
        analysis_added = False

        for child in child_runs:
            if child.workflow_stage != "Analysis":
                stages.append(
                    {
                        "id": child.workflow_stage
                        or f"step-{child.workflow_step_order}",
                        "label": self._workflow_stage_label(
                            child.workflow_stage or ""
                        ),
                        "state": self._workflow_child_state(child),
                        "child_run_id": child.id,
                        "application": self._workflow_child_application(child),
                        "status": self._workflow_child_status(child),
                    }
                )
                continue
            if analysis_added:
                continue

            analysis_added = True
            child_states = [
                self._workflow_child_state(analysis_child)
                for analysis_child in analysis_children
            ]
            if child_states and all(state == "complete" for state in child_states):
                state = "complete"
                status = "complete"
            elif "active" in child_states:
                state = "active"
                status = "active"
            else:
                state = "pending"
                status = "pending"
            selected_child = next(
                (
                    analysis_child
                    for analysis_child, child_state in zip(
                        analysis_children, child_states
                    )
                    if child_state != "complete"
                ),
                analysis_children[-1],
            )
            stages.append(
                {
                    "id": "Analysis",
                    "label": self._workflow_stage_label("Analysis"),
                    "state": state,
                    "status": status,
                    "child_run_id": selected_child.id,
                    "child_run_ids": [
                        analysis_child.id
                        for analysis_child in analysis_children
                    ],
                    "application": self._workflow_child_application(
                        selected_child
                    ),
                    "applications": [
                        self._workflow_child_application(analysis_child)
                        for analysis_child in analysis_children
                    ],
                }
            )
        return stages

    def _with_visualization_stage(self, stages):
        if any(stage.get("id") == "Visualization" for stage in stages):
            return stages

        analysis_stages = [stage for stage in stages if stage.get("id") == "Analysis"]
        analysis_complete = bool(analysis_stages) and all(
            stage.get("state") == "complete"
            or stage.get("status") in ("complete", "imported")
            for stage in analysis_stages
        )
        source_stage = next(
            (
                stage
                for stage in reversed(analysis_stages)
                if stage.get("child_run_id")
            ),
            None,
        )
        return [
            *stages,
            {
                "id": "Visualization",
                "label": "Visualization",
                "state": "active" if analysis_complete else "pending",
                "status": "ready" if analysis_complete else "pending",
                "local_only": True,
                "source_child_run_id": (
                    source_stage.get("child_run_id") if source_stage else None
                ),
            },
        ]

    def _workflow_child_state(self, child_run):
        if child_run.workflow_step_status == "complete" or self._workflow_child_has_finished_execution(child_run):
            return "complete"
        if child_run.workflow_step_status in ("submitted", "running"):
            return "active"
        if child_run.executions.exists():
            return "active"
        return "pending"

    def _workflow_child_status(self, child_run):
        if child_run.workflow_step_status == "complete" or self._workflow_child_has_finished_execution(child_run):
            return "complete"
        return child_run.workflow_step_status or "pending"

    def _workflow_child_has_finished_execution(self, child_run):
        request = self.context.get("request")
        return bool(
            request
            and child_run.executions.exists()
            and child_run.are_all_executions_finished(request)
        )

    def _workflow_child_application(self, child_run):
        return (
            child_run.workflow_application
            or child_run.module_application
            or child_run.utility_application
            or ""
        )

    def _workflow_progress(self, run_instance, stages):
        metadata = run_instance.workflow_metadata or {}
        workflow_state = metadata.get("workflow_state") or "not_started"
        active_child_run_id = metadata.get("active_child_run_id")
        active_stage = next(
            (
                stage for stage in stages
                if active_child_run_id and stage.get("child_run_id") == active_child_run_id
            ),
            None,
        )
        if active_stage is None:
            active_stage = next(
                (stage for stage in stages if stage.get("state") == "active"),
                None,
            )
        if active_stage and workflow_state == "not_started":
            workflow_state = "running"

        return {
            "workflow_state": workflow_state,
            "active_child_run_id": active_stage.get("child_run_id") if active_stage else active_child_run_id,
            "active_child_label": active_stage.get("label") if active_stage else "",
            "active_child_status": active_stage.get("status") if active_stage else "",
        }

    def _default_code(self, run_instance):
        return "\n".join(
            [
                f"# {run_instance.name}",
                "Calculation_Type = MODULE",
                "EPOLYSCAT_Application_Module = EPOLYSCAT_DMAT",
                "Input_File = target",
            ]
        )

    def get_is_tutorial(self, run_instance: models.Run):
        request = self.context["request"]
            
        tutorial_view = models.View.tutorial_view()
            
        return tutorial_view in run_instance.views.all() and django_user(request) != tutorial_view.owner
            
#    def get_status(self, instance: models.Run):
#        request = self.context["request"]
#        if not instance.executions.exists():
#            return "Unsubmitted"
#        else:
#            # get the last execution and return it's status
#            latest_execution: models.RemoteExecution = instance.latest_execution
#            # If not finished, try to get application specific status
#            if not latest_execution.is_airavata_experiment_finished(request):
#                application_status = latest_execution.get_application_specific_status(
#                    request
#                )
#                if application_status is not None:
#                    return application_status
#            return latest_execution.get_airavata_experiment_status(request)

    def get_status(self, run_instance: models.Run):
        request = self.context["request"]
        if not run_instance.executions.exists():
            return "UNSUBMITTED"
        else:
            # get the last execution and return it's status
            latest_execution: models.RemoteExecution = run_instance.latest_execution
            # If not finished, try to get application specific status
            if not latest_execution.is_airavata_experiment_finished(request):
                # experiment: ExperimentModel = request.airavata_client.getExperiment(
                #     request.authz_token, latest_execution.airavata_experiment_id
                # )

                # application_status = experiment_util.intermediate_output.get_intermediate_output_process_status(
                #     request, experiment, "bsr_prep.log"
                # )

                application_status = request.airavata_client.getExperimentStatus(
                    request.authz_token, latest_execution.airavata_experiment_id
                )

                if application_status is not None:
                    state = application_status.state
                    status = "CREATED" if state == 0 else (
                        "VALIDATED" if state == 1 else
                        "SCHEDULED"     if state == 2 else
                        "LAUNCHED"      if state == 3 else
                        "EXECUTING"     if state == 4 else
                        "CANCELING"     if state == 5 else
                        "CANCELED"      if state == 6 else
                        "COMPLETED"     if state == 7 else
                        "FAILED"
                    )

                    return status

            return latest_execution.get_airavata_experiment_status(request)



    def get_job_status(self, run_instance: models.Run):
        request = self.context["request"]
            
        if not run_instance.executions.exists():
            return "UNSUBMITTED"
        else:
            # get the last execution and return it's status
            latest_execution: models.RemoteExecution = run_instance.latest_execution

            try:
                job_statuses = request.airavata_client.getJobStatuses(
                    request.authz_token, latest_execution.airavata_experiment_id
                )
        
                job_statuses_list = list(job_statuses.values());
            
                if len(job_statuses_list) > 0:
                    # gets the most recent status
                    job_statuses_list.sort(key=lambda status: status.timeOfStateChange, reverse=True)
                    state = job_statuses_list[0].jobState

                    status = "SUBMITTED"    if state == 0 else (
                        "QUEUED"            if state == 1 else
                        "ACTIVE"            if state == 2 else
                        "COMPLETED"         if state == 3 else
                        "CANCELED"          if state == 4 else
                        "FAILED"            if state == 5 else
                        "SUSPENDED"         if state == 6 else
                        "UNKNOWN"           if state == 7 else
                        "NON_CRITICAL_FAIL"     if state == 8 else "UNKNOWN_"
                    )

                    return status
                else:
                    return "NO STATUS"
            except:
                return "---"

    def get_job_id(self, run_instance: models.Run):
        request = self.context["request"]

        if not run_instance.executions.exists():
            return None
        else:
            # get the last execution and return it's status
            latest_execution: models.RemoteExecution = run_instance.latest_execution
            return latest_execution.get_job_id(request)



    def get_resource(self, instance):
        request = self.context["request"]
        if not instance.executions.exists():
            return ""
        else:
            # get the last execution and return it's status
            latest_execution = instance.latest_execution
            return latest_execution.resource_name

    def get_resource_short(self, instance):
        request = self.context["request"]
        if not instance.executions.exists():
            return ""
        else:
            # get the last execution and return it's status
            latest_execution = instance.latest_execution
            return latest_execution.resource_name_short

    def get_can_resubmit(self, instance):
        request = self.context["request"]
        job_id = instance.get_most_recent_job_id(request)
        all_finished = instance.are_all_executions_finished(request)
        return job_id is not None and all_finished

    def get_cancelable(self, instance: models.Run):
        request = self.context["request"]
        return instance.is_cancelable(request)

    def _create_inpc_file(self, instance, directedit):
        request = self.context["request"]
        directedit_file = StringIO(directedit)
        data_product = user_storage.save(
            request,
            instance.filepath,
            file=directedit_file,
            name="inpc",
            content_type="text/plain",
        )
        instance.inpc_data_product_uri = data_product.productUri
        instance.save()

    def _create_inpc_file_input_table(self, instance, input_table):
        request = self.context["request"]
        input_table_file = StringIO(
            self._create_inpc_string_from_input_table(input_table)
        )
        data_product = user_storage.save(
            request,
            instance.filepath,
            file=input_table_file,
            name="inpc",
            content_type="text/plain",
        )
        instance.inpc_data_product_uri = data_product.productUri
        instance.input_table = json.dumps(input_table)
        instance.save()

    def _create_inpc_string_from_input_table(self, input_table):
        input_table_file = StringIO()
        input_table_file.write("# --- uRecX: machine-generated by uRecX ---")
        for pag in input_table["pages"]:
            for sec in pag["sections"]:
                names = [item["name"] for item in sec["lines"][0]["items"]]
                head = "\n" + (sec["category"] + ": ") + ",".join(names)
                for nlin in range(len(sec["lines"])):
                    s = "\n"
                    for item in sec["lines"][nlin]["items"]:
                        val = item["value"]
                        if val == "FLAG_ONLY":
                            val = ""
                        if val == "OBSOLETE":
                            val = ""
                        if val.find(",") != -1 or val.find(":") != -1:
                            # enclose in quotation if contains comma
                            val = "'" + val + "'"
                        s += val + ","
                    # Only write the header if there are values and this is the first line
                    if s.replace(",", "").strip() != "":
                        if nlin == 0:
                            input_table_file.write(head)
                    # If there are values, write them, leaving off the final trailing comma
                    if s.replace(",", "").strip():
                        input_table_file.write(s[:-1])
                input_table_file.write("\n\n")
        # Rewind to the begin of the file before trying to read it
        input_table_file.seek(0)
        return input_table_file.read()


class ExperimentSerializer(serializers.ModelSerializer):
    owner = serializers.SlugRelatedField(slug_field="username", read_only=True)
    run_count = serializers.SerializerMethodField()
    active_run_count = serializers.SerializerMethodField()
    description = serializers.CharField(allow_blank=True)
    name = serializers.CharField(
        required=True,
        validators=[UniqueToUserValidator(models.Experiment.objects.all(), "owner")],
    )

    class Meta:
        model = models.Experiment
        fields = (
            "id",
            "name",
            "description",
            "owner",
            "created",
            "updated",
            "deleted",
            "airavata_project_id",
            "run_count",
            "active_run_count",
        )
        read_only_fields = ("deleted", "airavata_project_id")

    def get_run_count(self, obj):
        return obj.runs.count()

    def get_active_run_count(self, obj):
        return obj.runs.filter(Q(deleted=False)).count()

    @transaction.atomic
    def create(self, validated_data):
        request = self.context["request"]
        experiment = models.Experiment.objects.create(
            **validated_data,
            owner=django_user(request),
        )
        try:
            experiment.create_airavata_project(request)
            experiment.save()
        except Exception:
            # Gateway connectivity is not required to create an experiment;
            # airavata_project_id stays null and can be provisioned later.
            logger.warning(
                "Failed to create Airavata project for experiment %s",
                experiment.name, exc_info=True,
            )
        return experiment

    @transaction.atomic
    def update(self, instance, validated_data):
        request = self.context["request"]
        # Don't allow updating name, since it must match the root name
        instance.description = validated_data["description"]
        instance.save()
        experiment = instance
        # For data migration, create an airavata project if there isn't one yet
        if experiment.airavata_project_id is None:
            try:
                experiment.create_airavata_project(request)
                experiment.save()
            except Exception:
                logger.warning(
                    "Failed to create Airavata project for experiment %s",
                    experiment.name, exc_info=True,
                )
        return experiment


class RunIdRelatedField(serializers.PrimaryKeyRelatedField):
    def get_queryset(self):
        request = self.context["request"]
        return models.Run.filter_by_user(request)

class PlotfileSerializer(serializers.Serializer):
    prefix = serializers.CharField(max_length=50, allow_blank=True)
    data_product_uri = serializers.CharField(max_length=200)


class PlotParametersIdRelatedField(serializers.PrimaryKeyRelatedField):
    def get_queryset(self):
        request = self.context["request"]
        return models.PlotParameters.filter_by_user(request)


class PlotParametersSerializer(serializers.ModelSerializer):
    xaxis = serializers.CharField(default="", allow_blank=True)
    yaxes = serializers.CharField(default="", allow_blank=True)
    flags = serializers.CharField(default="", allow_blank=True)

    class Meta:
        model = models.PlotParameters
        fields = (
            "id",
            "xaxis",
            "yaxes",
            "flags",
            "created",
            "last_use",
        )

    def create(self, validated_data):
        request = self.context["request"]
        plot_parameters, created = models.PlotParameters.objects.get_or_create(
            **validated_data,
            owner=django_user(request),
        )
        return plot_parameters

    def validate(self, attrs):
        attrs = super().validate(attrs)
        xaxis = attrs.get("xaxis", "")
        yaxes = attrs.get("yaxes", "")
        if xaxis and not yaxes:
            raise serializers.ValidationError(
                {"yaxes": ["yaxes is required when xaxis is also specified"]}
            )
        if yaxes and not xaxis:
            raise serializers.ValidationError(
                {"xaxis": ["xaxis is required when yaxes is also specified"]}
            )
        return attrs


class PlotSerializer(serializers.Serializer):
    runs = RunIdRelatedField(many=True)
    plotfile = serializers.CharField(max_length=20, required=False, allow_blank=True)
    plotfiles = PlotfileSerializer(many=True)
    plot_parameters = PlotParametersSerializer(required=False)
    plot_parameters_id = PlotParametersIdRelatedField(required=False)

    def validate(self, attrs):
        if "plot_parameters" not in attrs and "plot_parameters_id" not in attrs:
            raise serializers.ValidationError(
                "One of plot_parameters or plot_parameters_id is required"
            )
        return attrs


class ListInputsSerializer(serializers.Serializer):
    runs = RunIdRelatedField(many=True)


class DiffInputsSerializer(serializers.Serializer):
    runs = RunIdRelatedField(many=True)


class PlotablesSerializer(serializers.Serializer):
    runs = RunIdRelatedField(many=True)


class AddRemoveRunsSerializer(serializers.Serializer):
    runs = RunIdRelatedField(many=True)


class ViewSerializer(serializers.ModelSerializer):
    run_count = serializers.SerializerMethodField()
    active_run_count = serializers.SerializerMethodField()
    owner = serializers.SlugRelatedField(slug_field="username", read_only=True)
    runs = serializers.SerializerMethodField()

    class Meta:
        model = models.View
        fields = ("id", "name", "owner", "created", "updated", "deleted", "type", "run_count", "runs", "active_run_count", 'is_tutorial')
        #read_only_fields = ("owner", "created", "updated", "deleted", "type")

    def get_runs(self, view_instance: models.View):
        runs = filter(
            lambda run: any(map(lambda view: view==view_instance,run.views.all())),
            models.Run.objects.all()
        )

        return RunSerializer(runs, many=True, context={'request': self.context['request']}).data

    def get_run_count(self, view_instance: models.View):
        return len(self.get_runs(view_instance))

#//    def get_run_count(self, obj):
#//        return obj.runs.exclude(experiment__owner=None).count()

    def get_active_run_count(self, obj):
        return obj.runs.exclude(experiment__owner=None).filter(Q(deleted=False)).count()

    @transaction.atomic
    def create(self, validated_data):
        request = self.context["request"]
        view = models.View.objects.create(
            **validated_data,
            type="user-defined",
            owner=django_user(request),
        )
        view.save()
        return view

    @transaction.atomic
    def update(self, instance, validated_data):
        if instance.type == "user-defined":
            instance.name = validated_data.pop("name", instance.name)

        instance.save()
        return instance
