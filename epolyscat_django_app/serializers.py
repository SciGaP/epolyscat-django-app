import json
import os
from io import StringIO
from urllib.parse import urlencode

from airavata_django_portal_sdk import user_storage
from django.db import transaction
from django.db.models import Q
from django.utils.text import get_valid_filename
from rest_framework import reverse, serializers, validators

from epolyscat_django_app import models


class UniqueToUserValidator(validators.UniqueValidator):
    requires_context = True

    def __init__(self, queryset, user_field, message=None, lookup="exact"):
        self.user_field = user_field
        super().__init__(queryset, message=message, lookup=lookup)

    def __call__(self, value, serializer_field):
        self.user = serializer_field.context["request"].user
        return super().__call__(value, serializer_field)

    def filter_queryset(self, value, queryset, field_name):
        # filter by current user
        queryset = queryset.filter(**{self.user_field: self.user})
        return super().filter_queryset(value, queryset, field_name)


class RunSerializer(serializers.ModelSerializer):
    root = serializers.CharField(max_length=100, required=False)
    directedit = serializers.CharField(
        style={"base_template": "textarea.html"},
        allow_blank=True,
        write_only=True,
        required=False,
    )
    inpc_download_url = serializers.SerializerMethodField()
    status = serializers.SerializerMethodField()
    resource = serializers.SerializerMethodField()
    resource_short = serializers.SerializerMethodField()
    executions = serializers.SlugRelatedField(
        slug_field="airavata_experiment_id", read_only=True, many=True
    )
    input_table = serializers.JSONField(allow_null=True, required=False)
    can_resubmit = serializers.SerializerMethodField()
    cancelable = serializers.SerializerMethodField()

    class Meta:
        model = models.Run
        fields = (
            "id",
            "name",
            "number",
            "root",
            "experiment",
            "created",
            "updated",
            "deleted",
            "directedit",
            "inpc_download_url",
            "group_resource_profile_id",
            "compute_resource_id",
            "queue_name",
            "core_count",
            "node_count",
            "walltime_limit",
            "total_physical_memory",
            "status",
            "resource",
            "resource_short",
            "executions",
            "input_table",
            "can_resubmit",
            "cancelable",
        )
        read_only_fields = ("deleted", "number", "experiment", "name")

    def to_representation(self, instance):
        rep = super().to_representation(instance)
        rep["root"] = instance.root.root
        if instance.input_table is not None:
            rep["input_table"] = json.loads(instance.input_table)
        return rep

    @transaction.atomic
    def create(self, validated_data):
        request = self.context["request"]
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

    def get_status(self, instance: models.Run):
        request = self.context["request"]
        if not instance.executions.exists():
            return "Unsubmitted"
        else:
            # get the last execution and return it's status
            latest_execution: models.RemoteExecution = instance.latest_execution
            # If not finished, try to get application specific status
            if not latest_execution.is_airavata_experiment_finished(request):
                application_status = latest_execution.get_application_specific_status(
                    request
                )
                if application_status is not None:
                    return application_status
            return latest_execution.get_airavata_experiment_status(request)

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
    root = serializers.SlugRelatedField(slug_field="root", read_only=True)

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
            "root",
        )
        read_only_fields = ("deleted", "airavata_project_id", "root")

    def get_run_count(self, obj):
        return obj.runs.count()

    def get_active_run_count(self, obj):
        return obj.runs.filter(Q(deleted=False)).count()

    @transaction.atomic
    def create(self, validated_data):
        request = self.context["request"]
        name = get_valid_filename(validated_data.pop("name"))
        root = models.RunsRoot.objects.create(root=name, owner=request.user)
        experiment = models.Experiment.objects.create(
            **validated_data,
            owner=request.user,
            root=root,
            name=name,
        )
        experiment.create_airavata_project(request)
        experiment.save()
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
            experiment.create_airavata_project(request)
            experiment.save()
        return experiment


class RunIdRelatedField(serializers.PrimaryKeyRelatedField):
    def get_queryset(self):
        request = self.context["request"]
        return models.Run.filter_by_user(request)


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
            owner=request.user,
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
    plotfile = serializers.CharField(max_length=20)
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

    class Meta:
        model = models.View
        fields = (
            "id",
            "name",
            "owner",
            "created",
            "updated",
            "deleted",
            "type",
            "run_count",
            "active_run_count",
        )
        read_only_fields = ("owner", "created", "updated", "deleted", "type")

    def get_run_count(self, obj):
        return obj.runs.exclude(experiment__owner=None).count()

    def get_active_run_count(self, obj):
        return obj.runs.exclude(experiment__owner=None).filter(Q(deleted=False)).count()

    @transaction.atomic
    def create(self, validated_data):
        request = self.context["request"]
        view = models.View.objects.create(
            **validated_data,
            type="user-defined",
            owner=request.user,
        )
        view.save()
        return view

    @transaction.atomic
    def update(self, instance, validated_data):
        if instance.type == "user-defined":
            instance.name = validated_data.pop("name", instance.name)

        instance.save()
        return instance
