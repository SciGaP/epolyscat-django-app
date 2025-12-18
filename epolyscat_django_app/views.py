import base64
import io
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import typing
from contextlib import contextmanager
from copy import deepcopy
from pathlib import Path

from airavata.model.application.io.ttypes import DataType
from airavata.model.experiment.ttypes import ExperimentModel, UserConfigurationDataModel
from airavata.model.scheduling.ttypes import ComputationalResourceSchedulingModel
from airavata.model.status.ttypes import ExperimentState
from airavata_django_portal_sdk import experiment_util, user_storage
from django.apps import apps
from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.core.exceptions import ObjectDoesNotExist
from django.db import transaction
from django.db.models import Count, Max, Q, Sum
from django.db.models.functions import Coalesce
from django.http import FileResponse
from django.shortcuts import render
from django.utils import timezone
from rest_framework import (
    exceptions,
    pagination,
    permissions,
    response,
    status,
    viewsets,
)
from rest_framework.decorators import action, api_view
from rest_framework.response import Response

from epolyscat_django_app import models, serializers
from epolyscat_django_app.epolyscat_utils import Linp, is_empty

BASE_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
logger = logging.getLogger(__name__)

# Create your views here.


@login_required
def home(request):

    # Example code: Airavata API client
    # Make calls to the Airavata API from your view, for example:
    #
    # experiments = request.airavata_client.searchExperiments(
    #        request.authz_token, settings.GATEWAY_ID, request.user.username, filters={},
    #        limit=20, offset=0)
    #
    # The authorization token is always the first argument of Airavata API calls
    # and is available as 'request.authz_token'. Some API methods require a
    # 'gatewayID' argument and that is available on the Django settings object
    # as 'settings.GATEWAY_ID'.
    # For documentation on other Airavata API methods, see
    # https://docs.airavata.org/en/master/technical-documentation/airavata-api/.
    # The Airavata Django Portal uses the Airavata Python Client SDK:
    # https://github.com/apache/airavata/tree/master/airavata-api/airavata-client-sdks/airavata-python-sdk

    # Example code: user_storage module
    # In your Django views, you can make calls to the user_storage module to manage a user's files in the gateway
    #
    # user_storage.listdir(request, "")  # lists the user's home directory
    # user_storage.open_file(request, data_product_uri=...)  # open's a file for a given data_product_uri
    # user_storage.save(request, "path/in/user/storage", file)  # save a file to a path in the user's storage
    #
    # For more information as well as other user_storage functions, see https://airavata-django-portal-sdk.readthedocs.io/en/latest/

    return render(
        request,
        "epolyscat_django_app/application.html",
        {"project_name": "ePolyScat-Django-App"},
    )


class IsOwner(permissions.BasePermission):
    def has_object_permission(self, request, view, obj):
        # Only the owner has write access
        return request.user == obj.owner

class IsOwnerOrReadOnly(permissions.BasePermission):
    def has_object_permission(self, request, view, obj):
        if request.method in permissions.SAFE_METHODS:
            return True
        # Only the owner has write access
        return request.user == obj.owner



class Pagination(pagination.PageNumberPagination):
    page_size = 10
    page_size_query_param = "page_size"


class ExperimentViewSet(viewsets.ModelViewSet):
    serializer_class = serializers.ExperimentSerializer
    permission_classes = [permissions.IsAuthenticated, IsOwnerOrReadOnly]
    pagination_class = Pagination

    def get_queryset(self):
        request = self.request
        # user_projects = self.request.airavata_client.getUserProjects(
        #     request.authz_token, settings.GATEWAY_ID, request.user.username, -1, 0
        # )
        # project_ids = list(map(lambda p: p.projectID, user_projects))
        return (
            # Returns Experiments where user is Experiment owner or Experiment is shared via project
            models.Experiment.objects.filter(
                Q(deleted=False)
                & Q(owner=request.user)
                # & (Q(owner=request.user) | Q(airavata_project_id__in=project_ids))
            )
            .annotate(recent_activity_date=Coalesce(Max("runs__updated"), "updated"))
            .order_by("-recent_activity_date")
        )

    def perform_destroy(self, instance):
        instance.deleted = True
        instance.save()

    @action(detail=False, methods=["get"])
    def statistics(self, request):
        "Counts of undeleted experiments and runs"
        queryset = self.get_queryset()
        experiments_count = queryset.count()
        # Sum each experiments undeleted runs
        runs_count = queryset.annotate(
            num_runs=Count("runs", filter=Q(runs__deleted=False))
        ).aggregate(Sum("num_runs"))["num_runs__sum"]
        return Response(
            {
                "experiments_count": 0
                if experiments_count is None
                else experiments_count,
                "runs_count": 0 if runs_count is None else runs_count,
            }
        )


class RunViewSet(viewsets.ModelViewSet):
    serializer_class = serializers.RunSerializer
    permission_classes = [permissions.IsAuthenticated, IsOwnerOrReadOnly]
    pagination_class = Pagination

    DATA_TYPE_TO_FILENAME = {
        "STDOUT": "stdout",
        "STDERR": "stderr",
    }

    def get_queryset(self):
        request = self.request
        _queryset = models.Run.filter_by_user(request)
        _queryset = _queryset.filter(deleted=False)

        if self.request.query_params.get("experiment"):
            _queryset = _queryset.filter(
                Q(experiment=self.request.query_params.get("experiment"))
            )

        show_tutorial_runs = False
        if self.request.query_params.get("viewId"):
            view = models.View.filter_by_user(request).filter(
                id=self.request.query_params.get("viewId")
            )
            if len(view) > 0 and view[0].type == "unsubmitted":
                view[0].populate_unsubmitted_runs(request=request)
            if len(view) > 0 and view[0].type == "tutorial":
              show_tutorial_runs = True
              _queryset = _queryset.filter(
              Q(views=self.request.query_params.get("viewId"))
            )
        # Exclude tutorials from the general listing unless we're listing tutorials
        if self.action == "list" and not show_tutorial_runs:
            _queryset = _queryset.exclude(experiment__owner=None)

        return _queryset

        '''
        tutorial_view_id = models.View.tutorial_view().id
 
        return (
            # Returns Runs owned by the user
            models.Run.objects.filter(
                Q(owner=self.request.user) | Q(views__id__contains=tutorial_view_id)
            )
        )
        '''
    '''
        request = self.request
        _queryset = models.Run.filter_by_user(request)
        _queryset = _queryset.filter(deleted=False)

        if self.request.query_params.get("experiment"):
            _queryset = _queryset.filter(
                Q(experiment=self.request.query_params.get("experiment"))
            )

        show_tutorial_runs = False
        if self.request.query_params.get("viewId"):
            view = models.View.filter_by_user(request).filter(
                id=self.request.query_params.get("viewId")
            )
            if len(view) > 0 and view[0].type == "unsubmitted":
                view[0].populate_unsubmitted_runs(request=request)
            elif len(view) > 0 and view[0].type == "tutorial":
                show_tutorial_runs = True

            _queryset = _queryset.filter(
                Q(views=self.request.query_params.get("viewId"))
            )

        # Exclude tutorials from the general listing unless we're listing tutorials
        if self.action == "list" and not show_tutorial_runs:
            _queryset = _queryset.exclude(experiment__owner=None)

        #return _queryset
        return (
            # Returns Runs owned by the user
            models.Run.objects.filter(
                Q(owner=self.request.user) | Q(views__id__contains=tutorial_view_id)
            )
        )
    ''' 
    # Note: transaction.atomic doesn't affect user_storage calls, it is possible
    # that files are created/deleted and then the method later errors out
    # with the database changes being undone
    @transaction.atomic
    def create(self, request, *args, **kwargs):
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        run = serializer.save(owner=request.user)

        user_storage.create_user_dir(request, dir_names=run.directory.split("/"))

        for input in request.data["inputs_data"]:
            self._create_input(request, run, input)

        if "viewIds" in request.data:
            for viewId in request.data["viewIds"]:
                if int(viewId) < 0:
                    try:
                        tutorial_view = models.View.tutorial_view()
                        run.views.add(tutorial_view)
                    except models.View.DoesNotExist:
                        pass
                else:
                    view = models.View.objects.get(pk=int(viewId))
                    run.views.add(view)
   
        return Response(serializer.data)
    @transaction.atomic
    def update(self, request, *args, **kwargs):
        partial = kwargs.pop('partial', False)
        run = self.get_object()
        serializer = self.get_serializer(run, data=request.data, partial=partial)
        serializer.is_valid(raise_exception=True)
        self.perform_update(serializer)

        if run.owner != request.user:
            raise Exception("You can only update a run that you own")

        for updated_input in request.data["inputs_data"]:
            self._update_input(request, run, updated_input)

        return Response(serializer.data)


    def _create_input(self, request, run_instance, input):
        if input["type"] == "files":
            new_input = models.Input.objects.create(
                type="files",
                run=run_instance,
                name=input["name"],
                value=None
            )

            # Right now I'm assuming that URI_COLLECTIONS are supposed
            # to be saved with product uri's seperated by commas
            for file_data in input["files"]:
                self._save_file(request, run_instance, file_data, new_input)
        else:
            # Ensures that if a run type is previously defined as something else, it gets overidden when updated
            if input["name"] in ["EPOLYSCAT_Application_Module", "Application_Utility", "Application_Workflow"]:
                matching_inputs = list(filter(lambda input:
                    input.name in ["EPOLYSCAT_Application_Module", "Application_Utility", "Application_Workflow"],
                    run_instance.inputs.all()
                ));

                for matching_input in matching_inputs:
                    matching_input.delete()

            models.Input.objects.create(
                type=input["type"],
                run=run_instance,
                name=input["name"],
                value=input["value"]
            )

    def _update_input(self, request, run_instance, updated_input):
        matching_inputs = list(filter(lambda input:
            input.type == updated_input["type"] and
            input.name == updated_input["name"],
            run_instance.inputs.all()
        ))

        if not matching_inputs:
            self._create_input(request, run_instance, updated_input)
        else:
            old_input = matching_inputs[0]

            if updated_input["type"] == "files":
                for updated_file in updated_input["files"]:
                    matching_files = list(filter(lambda file:
                        file.name == updated_file["name"],
                        old_input.files.all()
                    ))

                    # file_exists = "dataProductURI" in updated_file and user_storage.exists(
                    #     request,
                    #     data_product_uri=updated_file["dataProductURI"]
                    # )

                    if not matching_files:
                        self._save_file(request, run_instance, updated_file, old_input)
                    else:
                        old_file = matching_files[0]

                        user_storage.delete(
                            request,
                            data_product_uri=old_file.data_product_uri
                        )

                        old_file.delete()

                        if not updated_file["deleted"]:
                            self._save_file(request, run_instance, updated_file, old_input)
            else:
                old_input.value = updated_input["value"]

            old_input.save()

    def _save_file(self, request, run_instance, file_data, input):
        if "deleted" not in file_data or not file_data["deleted"]:
            if "contents" not in file_data or file_data["contents"] == None:
                if not "dataProductURI" in file_data and "data-product-uri" in file_data:
                    file_data["dataProductURI"] = file_data["data-product-uri"]
                file = user_storage.open_file(
                    request, data_product_uri=file_data["dataProductURI"]
                )
                content_type = user_storage.get_data_product_metadata(
                    request, data_product_uri=file_data["dataProductURI"]
                )["mime_type"]
            elif file_data["isPlaintext"]:
                file = StringIO(file_data["contents"])
                content_type = "text/plain"
            else:
                file = BytesIO(base64.b64decode(file_data["contents"]))
                content_type = ""

            saved_file = user_storage.save(
                request,
                run_instance.directory,
                file,
                name=file_data["name"],
                content_type=content_type
            )

            models.File.objects.create(
                name=file_data["name"],
                data_product_uri=saved_file.productUri,
                input=input
            )


    @transaction.atomic
    def destroy(self, request, *args, **kwargs):
        run: models.Run = self.get_object()
        delete_associated = False if request.GET["deleteAssociated"] == "false" else True

        if run.owner != request.user:
            raise Exception("You can only delete a run that you own")

        if delete_associated and user_storage.dir_exists(request, run.directory):
            user_storage.delete_dir(request, run.directory)

        return super().destroy(request, *args, **kwargs)

    @transaction.atomic
    @action(detail=True, methods=["POST"])
    def clone(self, request, pk=None):
        run: models.Run = self.get_object()

        serializer = self.get_serializer(run)
        response = serializer.data

        del response["id"]
        del response["created"]
        del response["updated"]
        del response["deleted"]
        del response["executions"]

        response["name"] = f"{response['name']} CLONE"
        response["status"] = "UNSUBMITTED"
        response["job_status"] = "UNSUBMITTED"

        return Response(response)

    def perform_destroy(self, instance):
        instance.deleted = True
        instance.save()

    @action(detail=True, methods=["post"])
    def submit(self, request, pk=None):
        # Update the instance
        run: models.Run = self.get_object()
        serializer = self.get_serializer(run, data=request.data, partial=True)
        serializer.is_valid(raise_exception=True)
        serializer.save()

        if run.owner != request.user:
            raise Exception("You can only submit a run that you own")

        # change to api call
        app_interface_id = self._get_eployscat_app_interface_id(request)

        inputs = {};

        for input in run.inputs.all():
            if input.type == "files":
                data_product_uris = []

                for file in input.files.all():
                    opened_file = user_storage.open_file(
                        request, data_product_uri=file.data_product_uri
                    )
                    data_product_uris.append(user_storage.save_input_file(
                        request, opened_file, name=file.name
                    ).productUri)

                # still not sure if this is what URI_COLLECTION is expecting
                inputs[input.name] = ",".join(data_product_uris)
            else:
                inputs[input.name] = input.value

        self._create_remote_execution(request, run, inputs, app_interface_id, serializer.data["is_tutorial"])

        serializer = self.get_serializer(run)
        return Response(serializer.data)

    def status(self, request, pk=None):
        run: models.Run = self.get_object()

        return Response(run.status)


        # Copy the inpc file from the run directory into the experiment
        # directory by setting it as an input file
        inpc_file = user_storage.open_file(
            request, data_product_uri=run.inpc_data_product_uri
        )
        inpc_file_input_value = user_storage.save_input_file(
            request, inpc_file, name="inpc.inp"
        )

        input_values = {"Input-File": inpc_file_input_value.productUri}

        self._create_remote_execution(request, run, app_interface_id, input_values)

        serializer = self.get_serializer(run)
        return Response(serializer.data)

    @action(detail=True, methods=["post"])
    def resubmit(self, request, pk=None):
        run: models.Run = self.get_object()
        serializer = self.get_serializer(run, data=request.data)
        serializer.is_valid(raise_exception=True)
        serializer.save()

        # search through previous executions and find one with a job id
        job_id = run.get_most_recent_job_id(request)
        if job_id is None:
            raise Exception(
                f"Cannot resubmit run {run.id}, has no executions with a job"
            )

        app_interface_id = self._get_eployscat_app_interface_id(request)

        # prepare inputs to remote execution and call self._create_remote_execution
        inpc_file = user_storage.open_file(
            request, data_product_uri=run.inpc_data_product_uri
        )
        inpc_file_input_value = user_storage.save_input_file(
            request, inpc_file, name="inpc.inp"
        )
        input_values = {
            "Input-File": inpc_file_input_value.productUri,
            "Previous_JobID_Restart": job_id,
        }

        self._create_remote_execution(request, run, app_interface_id, input_values)

        serializer = self.get_serializer(run)
        return Response(serializer.data)

    def _create_remote_execution(
        self,
        request,
        run: models.Run,
        inputs: typing.Dict[str, str],
        app_interface_id: str,
        input_values: typing.Dict[str, str],
        is_tutorial: bool
    ) -> models.RemoteExecution:

        # Make sure that there aren't any currently running executions
        if not run.are_all_executions_finished(request):
            raise Exception("Run already has a currently running execution")

        # create experiment
        experiment = ExperimentModel()
        experiment.experimentName = f"{run.root.root}/{run.number} execution number {run.executions.count() + 1}"
        application_interface = request.airavata_client.getApplicationInterface(
            request.authz_token, app_interface_id
        )
        experiment.experimentInputs = application_interface.applicationInputs.copy()
        experiment.experimentOutputs = application_interface.applicationOutputs.copy()
        experiment.executionId = app_interface_id
        if run.experiment.airavata_project_id is None:
            run.experiment.create_airavata_project(request)
            run.experiment.save()
        experiment.projectId = run.experiment.airavata_project_id
        experiment.gatewayId = settings.GATEWAY_ID
        experiment.userName = request.user.username
        #ucd = UserConfigurationDataModel()
        #ucd.groupResourceProfileId = run.group_resource_profile_id
        #ucd.shareExperimentPublicly=is_tutorial

        #experiment.userConfigurationData = ucd
        experiment.userConfigurationData = UserConfigurationDataModel(
            groupResourceProfileId=run.group_resource_profile_id,
            shareExperimentPublicly=is_tutorial,
            computationalResourceScheduling=ComputationalResourceSchedulingModel(
                resourceHostId=run.compute_resource_id,
                totalCPUCount=run.core_count,
                nodeCount=run.node_count,
                wallTimeLimit=run.walltime_limit,
                queueName=run.queue_name
            )
        )
        #crs = ComputationalResourceSchedulingModel()
        #crs.resourceHostId = run.compute_resource_id
        #crs.totalCPUCount = run.core_count
        #crs.nodeCount = run.node_count
        #crs.wallTimeLimit = run.walltime_limit
        #crs.queueName = run.queue_name
        #experiment.userConfigurationData.computationalResourceScheduling = crs

        for inp in experiment.experimentInputs:
            if inp.name in input_values:
                inp.value = input_values[inp.name]

        # Save experiment
        experiment_id = request.airavata_client.createExperiment(
            request.authz_token, settings.GATEWAY_ID, experiment
        )
        # launch experiment
        experiment_util.launch(request, experiment_id)

        # add experiment to the run's executions
        compute_resource = request.airavata_client.getComputeResource(
            request.authz_token, run.compute_resource_id
        )
        return run.executions.create(
            airavata_experiment_id=experiment_id,
            resource_name=compute_resource.hostName,
        )
    @action(detail=True, methods=["PATCH"])
    def change_notification_settings(self, request, pk=None, *args, **kwargs):
        run = self.get_object()
        serializer = self.get_serializer(run, data=request.data, partial=True)
        serializer.is_valid(raise_exception=True)
        self.perform_update(serializer)

        execution = run.latest_execution

        experiment = request.airavata_client.getExperiment(
            request.authz_token, execution.airavata_experiment_id
        )
        if run.is_email_notification_on:
            experiment.emailAddresses = [request.user.email]
        else:
            experiment.emailAddresses = []

        request.airavata_client.updateExperiment(
            request.authz_token, execution.airavata_experiment_id, experiment
        )

        return Response(serializer.data)

    @action(detail=False, methods=["GET"])
    def tutorial_runs(self, request):
        tutorial_view = models.View.tutorial_view()
        runs = [run for run in models.Run.objects.all() if tutorial_view in run.views.all()]

        serializer = self.get_serializer(runs, many=True)

        return Response(serializer.data)


    def _get_eployscat_app_interface_id(self, request):
        app_module_id = getattr(settings, "EPOLYSCAT", {}).get(
            "EPOLYSCAT_APPLICATION_ID", "ePolyScat_940ab1c9-4ceb-431c-8595-c6246a195442"
        )
        all_app_interfaces = request.airavata_client.getAllApplicationInterfaces(
            request.authz_token, settings.GATEWAY_ID
        )
        app_interfaces = []
        for app_interface in all_app_interfaces:
            if not app_interface.applicationModules:
                continue
            if app_module_id in app_interface.applicationModules:
                app_interfaces.append(app_interface)
        if len(app_interfaces) == 1:
            app_interface_id = app_interfaces[0].applicationInterfaceId
        else:
            raise Exception(
                f"Could not figure out the applicationInterfaceId for app module {app_module_id}"
            )
        return app_interface_id

    @action(
        detail=True, methods=["post"], permission_classes=[permissions.IsAuthenticated]
    )
    def new(self, request, pk=None):
        run: models.Run = self.get_object()

        eployscat_settings = apps.get_app_config("epolyscat_django_app").APPLICATION_SETTINGS[
            "ePolyScat"
        ]

        # check if run has a linp file
        try:
            with open_run_file(request, run, "linp") as linp_file:
                linp = Linp(linp_file)
                include_values = True
        except FileNotFoundError:
            # otherwise use the master linp file
            linp = Linp(eployscat_settings["MASTER_LINP"])
            # since we're using the master linp file, we won't include any values from it
            include_values = False

        # load page/sections configuration for app
        INPUT_PAGES = deepcopy(eployscat_settings["INPUT_PAGES"])
        INPUT_PAGES["Other"] = deepcopy(eployscat_settings["ALL_INPUTS"])

        # Keep track of categories and names that have already been added
        all_sections = {}
        # generate JSON representation of page/sections with linp file
        # providing 'value', 'defaultValue' and 'documentation'
        input_table = {"pages": []}
        for input_page in INPUT_PAGES:
            sections = []
            page = {"name": input_page, "sections": sections}
            input_table["pages"].append(page)
            for input_section in INPUT_PAGES[input_page]:
                for input_category in input_section:
                    more_lines = True
                    line_num = 1
                    section = {"category": input_category, "lines": []}
                    section_names = []  # names to include in this section
                    # Figure out what category/names have already been included
                    for input_name in input_section[input_category]:
                        if input_category not in all_sections:
                            all_sections[input_category] = []
                        if input_name in all_sections[input_category]:
                            # Skip this input_name since it has been already included
                            continue
                        else:
                            all_sections[input_category].append(input_name)
                            section_names.append(input_name)
                    while more_lines and len(section_names) > 0:
                        items = []
                        for input_name in section_names:
                            linp_item = linp.item(
                                input_category, input_name, line=line_num
                            )
                            items.append(
                                {
                                    "name": input_name,
                                    "value": linp_item.inputValue
                                    if linp_item
                                    and include_values
                                    and not linp_item.empty
                                    else "",
                                    "defaultValue": linp_item.default
                                    if linp_item
                                    else None,
                                    "documentation": linp_item.docu
                                    if linp_item
                                    else None,
                                }
                            )
                        section["lines"].append({"items": items})
                        # Check to see if there are any more lines by
                        # checking if there are any items in the next line
                        current_line_num = line_num
                        for input_name in section_names:
                            if (
                                linp.item(input_category, input_name, line=line_num + 1)
                                is not None
                            ):
                                line_num = line_num + 1
                                logger.debug(
                                    f"For {input_category}:{input_name} incrementing line_num to {line_num}"
                                )
                                break
                        if not include_values or line_num == current_line_num:
                            more_lines = False
                    if len(section["lines"]) > 0:
                        sections.append(section)

        serializer = self.get_serializer(run)
        response = serializer.data
        # remove id, name, number, created, updated, deleted, executions
        del response["id"]
        del response["name"]
        del response["number"]
        del response["created"]
        del response["updated"]
        del response["deleted"]
        del response["executions"]
        # handle inpc_download_url for tutorials
        if run.is_tutorial:
            with open_run_file(request, run, "inpc") as inpc_file:
                inpc_dp = user_storage.save_input_file(
                    request, inpc_file, name="inpc", content_type="text/plain"
                )
                inpc_download_url = user_storage.get_download_url(
                    request, data_product=inpc_dp
                )
                response["inpc_download_url"] = inpc_download_url
        # otherwise, leave in inpc_download_url, this way frontend can retrieve
        # source run's input file

        # update status to Unsubmitted
        response["status"] = "Unsubmitted"

        # add input table
        response["input_table"] = input_table

        return Response(response)

    @action(methods=["get"], detail=True)
    def viewables(self, request, pk=None):
        run: models.Run = self.get_object()
        viewables = []

        epolyscat_settings = apps.get_app_config("epolyscat_django_app").APPLICATION_SETTINGS[
            "ePolyScat"
        ]
        for filename, description in epolyscat_settings["FILE_VIEWABLE"].items():
            if run_file_exists(request, run, filename):
                url = self.reverse_action(url_name="show-viewable", args=[pk, filename])
                viewables.append(dict(filename=filename, url=url))

        for data_type in self.DATA_TYPE_TO_FILENAME.keys():
            try:
                dp_uri = get_run_output_data_product_uri(
                    request, run, data_type=data_type
                )
                if dp_uri is not None:
                    filename = self.DATA_TYPE_TO_FILENAME[data_type]
                    url = self.reverse_action(
                        url_name="show-viewable", args=[pk, filename]
                    )
                    viewables.append(dict(filename=filename, url=url))
            except:
                logger.exception(f"Failed to check if {data_type} is available")

        return Response(viewables)

    @action(methods=["get"], detail=True, url_path="input-files")
    def input_files(self, request, pk=None):
        run: models.Run = self.get_object()
        input_files_list = []

        eployscat_settings = apps.get_app_config("epolyscat_django_app").APPLICATION_SETTINGS[
            "ePolyScat"
        ]
        for filename, description in eployscat_settings["FILE_INPUT"].items():
            if run_file_exists(request, run, filename):
                url = self.reverse_action(url_name="show-viewable", args=[pk, filename])
                input_files_list.append(dict(filename=filename, url=url))
        return Response(input_files_list)

    @action(methods=["get"], detail=True, url_path="viewables/(?P<filename>[\w]+)")
    def show_viewable(self, request, pk=None, filename: str = None):

        run: models.Run = self.get_object()
        try:
            # file may be experiment output file or a file from the run directory
            if filename in self.DATA_TYPE_TO_FILENAME.values():
                data_type = [
                    dt
                    for dt, fn in self.DATA_TYPE_TO_FILENAME.items()
                    if fn == filename
                ][0]
                dp_uri = get_run_output_data_product_uri(
                    request, run, data_type=data_type
                )
                f = user_storage.open_file(request, data_product_uri=dp_uri)
            else:
                f = open_run_file(request, run, filename)
            response = FileResponse(f, content_type="text/plain; charset=utf-8")
            return response
        except (ObjectDoesNotExist, FileNotFoundError) as e:
            raise exceptions.NotFound from e

    @action(methods=["put"], detail=True)
    def cancel(self, request, pk=None):
        run: models.Run = self.get_object()
        remote_execution: models.RemoteExecution = run.latest_execution
        if remote_execution is None:
            raise Exception("Run has no executions to cancel")
        remote_execution.cancel(request)
        serializer = self.get_serializer(run)
        return Response(serializer.data)


class PlotParametersViewSet(viewsets.ModelViewSet):
    serializer_class = serializers.PlotParametersSerializer
    permission_classes = [permissions.IsAuthenticated, IsOwner]

    def get_queryset(self):
        request = self.request
        return models.PlotParameters.filter_by_user(request).order_by("-last_use")


@api_view(["POST"])
def plot(request):
    "Returns dictionary with 'mime-type' and 'plot' as base64 encoded image."
    serializer = serializers.PlotSerializer(
        data=request.data, context={"request": request}
    )
    serializer.is_valid()
    serializer.is_valid(raise_exception=True)

    plot_command = [
        sys.executable,
        os.path.join(apps.get_app_config("epolyscat_django_app").SCRIPTS, "plot.py"),
    ]

    plotfiles = serializer.validated_data["plotfiles"]

    def map_plotfile_to_obj(plotfile):
        file = user_storage.open_file(request, data_product_uri=plotfile["data_product_uri"])
        sep = "" if plotfile['prefix'] == "" else "__"
        file_name = file.name.split("/")[-1]
        return { "file": file, "name": f"{plotfile['prefix']}{sep}{file_name}" }

    plotfile_objs = list(map(map_plotfile_to_obj, plotfiles))
    # with runs_dir(request, runs, [plotfile], ignore_missing=True) as tmp_runs_dir:
    with get_tmp_dir(request, plotfile_objs, False) as tmp_dir:
        plot_parameters = serializer.validated_data.get("plot_parameters", None)
        if plot_parameters is not None:
            obj, created = models.PlotParameters.objects.get_or_create(
                xaxis=plot_parameters["xaxis"],
                yaxes=plot_parameters["yaxes"],
                flags=plot_parameters["flags"],
                owner=request.user,
            )
            plot_parameters = obj
            if not created:
                plot_parameters.last_use = timezone.now()
                plot_parameters.save()
        else:
            plot_parameters = serializer.validated_data["plot_parameters_id"]
            plot_parameters.last_use = timezone.now()
            plot_parameters.save()

        xaxis = plot_parameters.xaxis
        yaxes = plot_parameters.yaxes
        flags = plot_parameters.flags
        cols = ""
        if xaxis:
            cols = xaxis + ":"
        if yaxes:
            cols += yaxes
        show_columns_command = plot_command.copy()
        for plotfile_obj in plotfile_objs:
            if cols:
                plot_command.append(
                    os.path.join(plotfile_obj["name"] + "[" + cols + "]")
                )
            else:
                plot_command.append(os.path.join(plotfile_obj["name"]))

        # For every run directory where the plotfile exists, add to the command
        # for run in runs:
        #     # command is executed in tmp_runs_dir and command line paths are
        #     # relative to that base directory
        #     rundir = os.path.join(run.root.root, run.number)
        #     if os.path.exists(os.path.join(tmp_runs_dir, rundir, plotfile)):
        #         if cols:
        #             plot_command.append(
        #                 os.path.join(rundir, plotfile + "[" + cols + "]")
        #             )
        #         else:
        #             plot_command.append(os.path.join(rundir, plotfile))
        #         # Don't include columns for showColumns in case they are out of range
        #         show_columns_command.append(os.path.join(rundir, plotfile))

        if flags:
            plot_command += flags.split()
        graph = "plot.png"
        plot_command.append("-plotfile=" + graph)
        plot_command.append("-batch")

        try:
            logger.debug(f"Running {plot_command}")
            process = subprocess.Popen(
                plot_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                cwd=tmp_dir,
            )
            returncode = process.wait()
            process_output = process.stdout.read().decode()
            logger.debug(f"plot.py returncode={returncode}")
            logger.debug(f"plot.py output={process_output}")

            data = {"mime-type": "image/png"}
            data["output"] = process_output
            if os.path.exists(os.path.join(tmp_dir, graph)):
                with open(os.path.join(tmp_dir, graph), "rb") as plotpng:
                    data["plot"] = base64.b64encode(plotpng.read()).decode("utf-8")

            # If there was an error, return the showColumns output
            if returncode != 0:
                show_columns_command.append("-batch")
                show_columns_command.append("-showColumns")
                logger.debug(f"Running showColumns {show_columns_command}")
                process = subprocess.Popen(
                    show_columns_command,
                    stdout=subprocess.PIPE,
                    # stderr=subprocess.STDOUT,
                    cwd=tmp_dir,
                )
                process.wait()
                process_output = process.stdout.read().decode()
                logger.debug(f"plot.py -showColumns returncode={process.returncode}")
                logger.debug(f"plot.py -showColumns output={process_output}")
                data["user_guidance"] = process_output

            return response.Response(data)
        except Exception as e:
            plot_parameters.delete()
            logger.exception(f"Failed to generate plot for {plot_command}")
            raise  # re-raise exception


@contextmanager
def get_tmp_dir(request, plotfile_objs, ignore_missing=False):
    "Copy run files into temporary directory, keeping the existing structure"
    with tempfile.TemporaryDirectory() as tmpdir:
        for plotfile in plotfile_objs:
            file = plotfile["file"]
            filename = plotfile["name"]

            if os.path.exists(os.path.join(tmpdir, filename)):
                i = 2

                while os.path.exists(os.path.join(tmpdir, filename + "__" + str(i))):
                    i += 1

                filename += "__" + str(i)

            with open(os.path.join(tmpdir, filename), "xb") as fcopy:
                shutil.copyfileobj(file, fcopy)

        yield tmpdir


    # with tempfile.TemporaryDirectory() as tmpdir:
    #     for run in runs:
    #         for filename in filenames:
    #             try:
    #                 with open_run_file(request, run, filename) as f:
    #                     rundir = os.path.join(tmpdir, run.root.root, run.number)
    #                     os.makedirs(rundir, exist_ok=True)
    #                     logger.debug(f"copying {run.filepath}/{filename} to {rundir}")
    #                     with open(os.path.join(rundir, filename), "wb") as fcopy:
    #                         shutil.copyfileobj(f, fcopy)
    #             except FileNotFoundError:
    #                 logger.debug(f"could not find {run.filepath}/{filename}")
    #                 if not ignore_missing:
    #                     raise
    #     yield tmpdir


def open_run_file(request, run: models.Run, filename: str):
    # Handle tutorial runs specially: all of their files are available within the app
    if False and run.is_tutorial:
        return open(BASE_DIR / run.filepath / filename, "rb")
    else:
        data_product_uri = user_run_file_exists(request, run, filename)
        if data_product_uri is not None:
            return user_storage.open_file(request, data_product_uri=data_product_uri)
        else:
            raise FileNotFoundError(f"{filename} does not exist in run {run.id}")


def run_file_exists(request, run: models.Run, filename: str) -> bool:
    if False and run.is_tutorial:
        return (BASE_DIR / run.directory / filename).exists()
    else:
        return user_run_file_exists(request, run, filename) is not None

def user_run_file_exists(request, run, filename):
    """Return data product uri for run file if it exists, else None."""

    # check to see if the file is already in the run directory, for backwards compatibility
    if run.owner == request.user:
        data_product_uri = user_storage.user_file_exists(
            request, os.path.join(run.directory, filename)
        )
        if data_product_uri is not None:
            logger.debug(f"Found {filename} in {run.filepath}")
            return data_product_uri

    # Find the most recent completed execution (Airavata experiment) if exists
    experiment_model = None
    for execution in run.executions.order_by("-created"):
        status_name = execution.get_airavata_experiment_status(request)
        experiment_state = ExperimentState._NAMES_TO_VALUES[status_name]
        if experiment_state == ExperimentState.COMPLETED:
            logger.debug(f"getExperiment({execution.airavata_experiment_id})")
            experiment_model = request.airavata_client.getExperiment(
                request.authz_token, execution.airavata_experiment_id
            )
            break
    if experiment_model is None:
        return None

    # Load the Modl_RunID file to find location of files
    # TODO: cache this information
    modl_runid_output = None
    for output in experiment_model.experimentOutputs:
        if output.name == "Modl_RunID":
            modl_runid_output = output
    if modl_runid_output is None or not user_storage.exists(
        request, data_product_uri=modl_runid_output.value
    ):
        raise Exception("Modl_RunID file is missing")
    modl_runid_file = user_storage.open_file(
        request, data_product_uri=modl_runid_output.value
    )
    model_runid = modl_runid_file.read().decode()
    m = re.match(r"(\S+) (\S+)", model_runid)
    if m is None:
        raise Exception(f"Invalid Modl_RunID file contents: {model_runid}")
    model, run_id = m.group(1, 2)


    # Check for the file in ARCHIVE/model/run_id/ directory
    data_product_uri = user_storage.user_file_exists(
        request,
        os.path.join("ARCHIVE", model, run_id, filename),
        experiment_id=experiment_model.experimentId,
    )
    return data_product_uri


@api_view(["GET"])
def api_settings(request):
    app_module_id = getattr(settings, "EPOLYSCAT", {}).get(
        "EPOLYSCAT_APPLICATION_ID",
        # "BSR:_B-Spline_atomic_R-matrix_code_9ae142cb-689f-4440-8d2d-e131f2891005"
        #"BSR3_82b15174-04a1-471e-82a3-33c77c8c6281"
        "ePolyScat_940ab1c9-4ceb-431c-8595-c6246a195442"
    )
    return response.Response({"EPOLYSCAT": {"EPOLYSCAT_APPLICATION_ID": app_module_id}})


def get_run_output_data_product_uri(request, run: models.Run, data_type: str):
    # Find most recent execution
    if run.executions.exists():
        most_recent_execution: models.RemoteExecution = run.executions.order_by(
            "-created"
        )[0]
    else:
        return None

    experiment_model: ExperimentModel = request.airavata_client.getExperiment(
        request.authz_token, most_recent_execution.airavata_experiment_id
    )
    # Find the output by data type
    output = None
    for output in experiment_model.experimentOutputs:
        output_type_name = DataType._VALUES_TO_NAMES[output.type]
        if output_type_name == data_type:
            output = output
            break
    # If experiment is finished, see if there is an experimentOutput available
    if most_recent_execution.is_airavata_experiment_finished(request):
        if output is not None or user_storage.exists(
            request, data_product_uri=output.value
        ):
            return output.value
        else:
            return None
    # otherwise, see if there is an intermediate output available
    else:
        data_products = (
            experiment_util.intermediate_output.get_intermediate_output_data_products(
                request, experiment_model, output.name
            )
        )
        if len(data_products) == 1 and user_storage.exists(
            request, data_product=data_products[0]
        ):
            return data_products[0].productUri
        else:
            return None


class ViewsViewSet(viewsets.ModelViewSet):
    serializer_class = serializers.ViewSerializer
    permission_classes = [permissions.IsAuthenticated, IsOwnerOrReadOnly]
    pagination_class = Pagination

    def get_queryset(self):
        #request = self.request
        return (
            # Returns Runs owned by the user
            models.View.objects.filter(
                Q(owner=self.request.user)
            )
        )

#        #queryset = models.View.filter_by_user(request)
#
#        if self.action == "list":
#            queryset = queryset.exclude(type="tutorial")
#
#        if self.action == "list" and "type" in self.request.query_params:
#            types = self.request.query_params.getlist("type")
#            logger.debug(f"filtering by types={types}")
#            queryset = queryset.filter(type__in=types)
#
#        queryset = (
#            queryset.filter(deleted=False)
#            .annotate(recent_activity_date=Coalesce(Max("runs__updated"), "updated"))
#            .order_by("-order", "-recent_activity_date")
#        )
#
#        return queryset

    @transaction.atomic
    def create(self, request, *args, **kwargs):
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        view = serializer.save(owner=request.user)

        for run in models.Run.objects.all():
            if run.id in request.data["runIds"]:
                run.views.add(view)

        return Response(serializer.data)

    @transaction.atomic
    def update(self, request, *args, **kwargs):
        partial = kwargs.pop('partial', False)
        view = self.get_object()
        serializer = self.get_serializer(view, data=request.data, partial=partial)
        serializer.is_valid(raise_exception=True)

        if request.data['overide']:
            for run in view.runs.all():
                run.views.remove(view)

        for run in models.Run.objects.all():
            if run.id in request.data["runIds"]:
                run.views.add(view)
        
                # if serializer.data["is_tutorial"]:
                    # for execution in run.executions:
                    #     experiment = request.airavata_client.getExperiment(
                    #         request.authz_token, execution.experiment_id
                    #     )

                    #     experiment.userConfigurationData.shareExperimentPublicly = True

        return Response(serializer.data)

''' 
     ###  changing tRecX to BSR model 

    def list(self, request, *args, **kwargs):
        models.View.create_default_views(request)
        if self.get_queryset().filter(type="unsubmitted").exists():
            unsubmitted_view = self.get_queryset().get(type="unsubmitted")
            unsubmitted_view.populate_unsubmitted_runs(request=request)
        return super().list(request, *args, **kwargs)

    def get_object(self):
        obj = super().get_object()
        if obj.type == "unsubmitted":
            obj.populate_unsubmitted_runs(request=self.request)
        return obj

    def perform_destroy(self, instance):
        if instance.type != "user-defined":
            raise Exception("Can only delete 'user-defined' views.")
        instance.deleted = True
        instance.save()

    @action(
        detail=True,
        methods=["put"],
        url_path="add-runs",
        serializer_class=serializers.AddRemoveRunsSerializer,
    )
    def add_runs(self, request, pk=None):
        view: models.View = self.get_object()
        if view.type not in ["default", "user-defined"]:
            raise Exception(
                "Can only add/remove runs from 'user-defined' and 'default' views."
            )

        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        view.runs.add(*serializer.validated_data["runs"])

        return Response(data=serializers.ViewSerializer(view).data)

    @action(
        detail=True,
        methods=["put"],
        url_path="remove-runs",
        serializer_class=serializers.AddRemoveRunsSerializer,
    )
    def remove_runs(self, request, pk=None):
        view: models.View = self.get_object()
        if view.type not in ["default", "user-defined"]:
            raise Exception(
                "Can only add/remove runs from 'user-defined' and 'default' views."
            )

        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        view.runs.remove(*serializer.validated_data["runs"])

        return Response(data=serializers.ViewSerializer(view).data)

    @action(methods=["GET"], detail=False)
    def tutorials(self, request):
        try:
            tutorials_view = models.View.objects.get(type="tutorial", owner=None)
            tutorials_view = self._refresh_tutorials_view(tutorials_view)
        except models.View.DoesNotExist:
            tutorials_view = self._create_tutorials_view()
        serializer = self.get_serializer(tutorials_view)
        return Response(serializer.data)

#    @transaction.atomic
#    def update(self, request, *args, **kwargs):
#        partial = kwargs.pop('partial', False)
#        view = self.get_object()
#        serializer = self.get_serializer(view, data=request.data, partial=partial)
#        serializer.is_valid(raise_exception=True)
#
#        if request.data['overide']:
#            for run in view.runs.all():
#                run.views.remove(view)
#
#        for run in models.Run.objects.all():
#            if run.id in request.data["runIds"]:
#                run.views.add(view)
#
#                # if serializer.data["is_tutorial"]:
#                    # for execution in run.executions:
#                    #     experiment = request.airavata_client.getExperiment(
#                    #         request.authz_token, execution.experiment_id
#                    #     )
#
#                    #     experiment.userConfigurationData.shareExperimentPublicly = True
#
#        return Response(serializer.data)

    @transaction.atomic
    def update(self, request, *args, **kwargs):
        partial = kwargs.pop('partial', False)
        run = self.get_object()
        serializer = self.get_serializer(run, data=request.data, partial=partial)
        serializer.is_valid(raise_exception=True)
        self.perform_update(serializer)

        if run.owner != request.user:
            raise Exception("You can only update a run that you own")

        for updated_input in request.data["inputs_data"]:
            self._update_input(request, run, updated_input)

        return Response(serializer.data)

    def _create_input(self, request, run_instance, input):
        if input["type"] == "files":
            new_input = models.Input.objects.create(
                type="files",
                run=run_instance,
                name=input["name"],
                value=None
            )

            # Right now I'm assuming that URI_COLLECTIONS are supposed
            # to be saved with product uri's seperated by commas
            for file_data in input["files"]:
                self._save_file(request, run_instance, file_data, new_input)
        else:
            # Ensures that if a run type is previously defined as something else, it gets overidden when updated
            if input["name"] in ["EPOLYSCAT_Application_Module", "Application_Utility", "Application_Workflow"]:
                matching_inputs = list(filter(lambda input:
                    input.name in ["EPOLYSCAT_Application_Module", "Application_Utility", "Application_Workflow"],
                    run_instance.inputs.all()
                ));
                for matching_input in matching_inputs:
                    matching_input.delete()

            models.Input.objects.create(
                type=input["type"],
                run=run_instance,
                name=input["name"],
                value=input["value"]
            )

    def _update_input(self, request, run_instance, updated_input):
        matching_inputs = list(filter(lambda input:
            input.type == updated_input["type"] and
            input.name == updated_input["name"],
            run_instance.inputs.all()
        ))

        if not matching_inputs:
            self._create_input(request, run_instance, updated_input)
        else:
            old_input = matching_inputs[0]

            if updated_input["type"] == "files":
                for updated_file in updated_input["files"]:
                    matching_files = list(filter(lambda file:
                        file.name == updated_file["name"],
                        old_input.files.all()
                    ))

                    # file_exists = "dataProductURI" in updated_file and user_storage.exists(
                    #     request,
                    #     data_product_uri=updated_file["dataProductURI"]
                    # )

                    if not matching_files:
                        self._save_file(request, run_instance, updated_file, old_input)
                    else:
                        old_file = matching_files[0]

                        user_storage.delete(
                            request,
                            data_product_uri=old_file.data_product_uri
                        )

                        old_file.delete()

                        if not updated_file["deleted"]:
                            self._save_file(request, run_instance, updated_file, old_input)
            else:
                old_input.value = updated_input["value"]

            old_input.save()


    def _save_file(self, request, run_instance, file_data, input):
        if "deleted" not in file_data or not file_data["deleted"]:
            if "contents" not in file_data or file_data["contents"] == None:
                if not "dataProductURI" in file_data and "data-product-uri" in file_data:
                    file_data["dataProductURI"] = file_data["data-product-uri"]
                file = user_storage.open_file(
                    request, data_product_uri=file_data["dataProductURI"]
                )
                content_type = user_storage.get_data_product_metadata(
                    request, data_product_uri=file_data["dataProductURI"]
                )["mime_type"]
            elif file_data["isPlaintext"]:
                file = StringIO(file_data["contents"])
                content_type = "text/plain"
            else:
                file = BytesIO(base64.b64decode(file_data["contents"]))
                content_type = ""

            saved_file = user_storage.save(
                request,
                run_instance.directory,
                file,
                name=file_data["name"],
                content_type=content_type
            )

            models.File.objects.create(
                name=file_data["name"],
                data_product_uri=saved_file.productUri,
                input=input
            )

    @transaction.atomic
    def destroy(self, request, *args, **kwargs):
        run: models.Run = self.get_object()
        delete_associated = False if request.GET["deleteAssociated"] == "false" else True

        if run.owner != request.user:
            raise Exception("You can only delete a run that you own")

        if delete_associated and user_storage.dir_exists(request, run.directory):
            user_storage.delete_dir(request, run.directory)

        return super().destroy(request, *args, **kwargs)


    @transaction.atomic
    @action(detail=True, methods=["POST"])
    def clone(self, request, pk=None):
        run: models.Run = self.get_object()

        serializer = self.get_serializer(run)
        response = serializer.data

        del response["id"]
        del response["created"]
        del response["updated"]
        del response["deleted"]
        del response["executions"]

        response["name"] = f"{response['name']} CLONE"
        response["status"] = "UNSUBMITTED"
        response["job_status"] = "UNSUBMITTED"

        return Response(response)

    @action(detail=True, methods=["GET"])
    def get_output_files(self, request, pk=None):
        lst = list(get_list_or_404(self.get_queryset(), pk=pk))
        run = self.get_object()
        output_files = []

        if len(run.executions.all()) > 0:
            most_recent_execution = run.executions.order_by(
                "-created"
            )[0]

            try:
                output_files = user_storage.list_experiment_dir(
                    request, most_recent_execution.airavata_experiment_id, path="ARCHIVE"
                )[1]
            except Exception:
                pass

        return Response(output_files)





    @transaction.atomic
    def _create_tutorials_view(self):

        logger.info("Creating Tutorials Experiment, Runs and View")
        tutorials_root = models.RunsRoot.objects.create(root="tutorial", owner=None)
        tutorials_experiment = models.Experiment.objects.create(
            name="Tutorials", owner=None, root=tutorials_root
        )
        tutorials_dir = Path(BASE_DIR, "data", "epolyscat", "tutorials")
        for tutorial_dir in tutorials_dir.iterdir():
            models.Run.objects.create(
                number=tutorial_dir.name,
                root=tutorials_root,
                experiment=tutorials_experiment,
                filepath=tutorial_dir.relative_to(BASE_DIR),
            )
        # Assume that the tutorials View also doesn't exist
        tutorials_view = models.View.objects.create(
            name="Tutorials", type="tutorial", owner=None
        )
        tutorials_view.runs.set(tutorials_experiment.runs.all())
        return tutorials_view

    @transaction.atomic
    def _refresh_tutorials_view(self, tutorials_view: models.View):
        # Figure out the tutorials experiment and root from one of the tutorial runs
        run: models.Run = tutorials_view.runs.first()
        tutorials_experiment: models.Experiment = run.experiment
        tutorials_root: models.RunsRoot = tutorials_experiment.root

        tutorials_dir = Path(BASE_DIR, "data", "epolyscat", "tutorials")
        for tutorial_dir in tutorials_dir.iterdir():
            run, created = models.Run.objects.get_or_create(
                number=tutorial_dir.name,
                root=tutorials_root,
                experiment=tutorials_experiment,
                defaults={
                    "filepath": tutorial_dir.relative_to(BASE_DIR),
                },
            )
            if created:
                logger.debug("Adding run %s to tutorials", str(run))
        tutorials_view.runs.set(tutorials_experiment.runs.all())
        return tutorials_view


@api_view(["POST"])
def plot(request):
    "Returns dictionary with 'mime-type' and 'plot' as base64 encoded image."
    serializer = serializers.PlotSerializer(
        data=request.data, context={"request": request}
    )
    serializer.is_valid(raise_exception=True)

    plot_command = [
        sys.executable,
        os.path.join(apps.get_app_config("epolyscat_django_app").SCRIPTS, "plot.py"),
    ]

    plotfile = serializer.validated_data["plotfile"]
    runs = serializer.validated_data["runs"]
    with runs_dir(request, runs, [plotfile], ignore_missing=True) as tmp_runs_dir:

        plot_parameters = serializer.validated_data.get("plot_parameters", None)
        if plot_parameters is not None:
            obj, created = models.PlotParameters.objects.get_or_create(
                xaxis=plot_parameters["xaxis"],
                yaxes=plot_parameters["yaxes"],
                flags=plot_parameters["flags"],
                owner=request.user,
            )
            plot_parameters = obj
            if not created:
                plot_parameters.last_use = timezone.now()
                plot_parameters.save()
        else:
            plot_parameters = serializer.validated_data["plot_parameters_id"]
            plot_parameters.last_use = timezone.now()
            plot_parameters.save()

        xaxis = plot_parameters.xaxis
        yaxes = plot_parameters.yaxes
        flags = plot_parameters.flags
        cols = ""
        if xaxis:
            cols = xaxis + ":"
        if yaxes:
            cols += yaxes

        show_columns_command = plot_command.copy()
        # For every run directory where the plotfile exists, add to the command
        for run in runs:
            # command is executed in tmp_runs_dir and command line paths are
            # relative to that base directory
            rundir = os.path.join(run.root.root, run.number)
            if os.path.exists(os.path.join(tmp_runs_dir, rundir, plotfile)):
                if cols:
                    plot_command.append(
                        os.path.join(rundir, plotfile + "[" + cols + "]")
                    )
                else:
                    plot_command.append(os.path.join(rundir, plotfile))
                # Don't include columns for showColumns in case they are out of range
                show_columns_command.append(os.path.join(rundir, plotfile))

        if flags:
            plot_command += flags.split()
        graph = "plot.png"
        plot_command.append("-plotfile=" + graph)
        plot_command.append("-batch")

        try:
            logger.debug(f"Running {plot_command}")
            process = subprocess.Popen(
                plot_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                cwd=tmp_runs_dir,
            )
            returncode = process.wait()
            process_output = process.stdout.read().decode()
            logger.debug(f"plot.py returncode={returncode}")
            logger.debug(f"plot.py output={process_output}")

            data = {"mime-type": "image/png"}
            data["output"] = process_output
            if os.path.exists(os.path.join(tmp_runs_dir, graph)):
                with open(os.path.join(tmp_runs_dir, graph), "rb") as plotpng:
                    data["plot"] = base64.b64encode(plotpng.read()).decode("utf-8")

            # If there was an error, return the showColumns output
            if returncode != 0:
                show_columns_command.append("-batch")
                show_columns_command.append("-showColumns")
                logger.debug(f"Running showColumns {show_columns_command}")
                process = subprocess.Popen(
                    show_columns_command,
                    stdout=subprocess.PIPE,
                    # stderr=subprocess.STDOUT,
                    cwd=tmp_runs_dir,
                )
                process.wait()
                process_output = process.stdout.read().decode()
                logger.debug(f"plot.py -showColumns returncode={process.returncode}")
                logger.debug(f"plot.py -showColumns output={process_output}")
                data["user_guidance"] = process_output

            return response.Response(data)
        except Exception as e:
            logger.exception(f"Failed to generate plot for {plot_command}")
            raise  # re-raise exception


@api_view(["POST"])
def list_inputs(request):
    "Returns dictionary with 'output' and 'stderr' as returned by lRuns"
    serializer = serializers.ListInputsSerializer(
        data=request.data, context={"request": request}
    )
    serializer.is_valid(raise_exception=True)

    runs = serializer.validated_data["runs"]

    try:
        result = _execute_lruns(request, runs)
        data = {"output": result.stdout, "stderr": result.stderr}
        if result.success:
            http_status = status.HTTP_200_OK
        else:
            http_status = status.HTTP_500_INTERNAL_SERVER_ERROR
        return response.Response(data, status=http_status)
    except FileNotFoundError as e:
        raise exceptions.NotFound(detail=str(e))


@api_view(["POST"])
def diff_inputs(request):
    "Returns dictionary with 'output' and 'stderr' as returned by lRuns"
    serializer = serializers.DiffInputsSerializer(
        data=request.data, context={"request": request}
    )
    serializer.is_valid(raise_exception=True)

    runs = serializer.validated_data["runs"]

    try:
        result = _execute_lruns(request, runs, diff=True)
        data = {"output": result.stdout, "stderr": result.stderr}
        if result.success:
            http_status = status.HTTP_200_OK
        else:
            http_status = status.HTTP_500_INTERNAL_SERVER_ERROR
        return response.Response(data, status=http_status)
    except FileNotFoundError as e:
        raise exceptions.NotFound(detail=str(e))


@api_view(["POST"])
def plotables(request):
    "Returns dictionary with 'filename' key and list of plotable files as value"

    serializer = serializers.PlotablesSerializer(
        data=request.data, context={"request": request}
    )
    serializer.is_valid(raise_exception=True)

    runs = serializer.validated_data["runs"]
    plotable_files = set()
    epolyscat_settings = apps.get_app_config("epolyscat_django_app").APPLICATION_SETTINGS[
        "ePolyScat"
    ]
    for filename, description in epolyscat_settings["FILE_PLOTABLE"].items():
        # Return filename if at least one run has the file
        for run in runs:
            try:
                with open_run_file(request, run, filename) as f:
                    if not is_empty(f):
                        plotable_files.add(filename)
            except FileNotFoundError:
                continue
    return response.Response({"filenames": plotable_files})


@api_view(["GET"])
def api_settings(request):
    app_module_id = getattr(settings, "EPOLYSCAT", {}).get(
        "EPOLYSCAT_APPLICATION_ID", "ePolyScat_940ab1c9-4ceb-431c-8595-c6246a195442"
    )
    return response.Response({"EPOLYSCAT": {"EPOLYSCAT_APPLICATION_ID": app_module_id}})
'''


class LRunsResult(typing.NamedTuple):
    stdout: str
    stderr: str
    success: bool


def _execute_lruns(request, runs, diff=False) -> LRunsResult:
    """Runs lRuns.py and returns stdout and stderr."""

    lruns_command = [
        sys.executable,
        os.path.join(apps.get_app_config("epolyscat_django_app").SCRIPTS, "lRuns.py"),
    ]

    with runs_dir(request, runs, ["inpc", "outf", "linp"]) as tmp_runs_dir:

        lruns_command.append(os.path.join(tmp_runs_dir, runs[0].root.root))
        lruns_command.append(",".join([run.number for run in runs]))
        if diff:
            lruns_command.append("-diff")

        # Run once to generate linp-extract
        with open(os.path.join(tmp_runs_dir, "lruns_out"), "w") as f:
            logger.debug(f"First time, to generate linp-extract: {lruns_command}")
            process = subprocess.Popen(lruns_command, stdout=f)
        returncode = process.wait()
        logger.debug(f"returncode={returncode}")

        # Run a second time to actually generate output
        with open(os.path.join(tmp_runs_dir, "lruns_out"), "w") as out, open(
            os.path.join(tmp_runs_dir, "lruns_err"), "w"
        ) as err:
            logger.debug(f"Second time: {lruns_command}")
            process = subprocess.Popen(lruns_command, stdout=out, stderr=err)
        returncode = process.wait()
        logger.debug(f"second returncode={returncode}")
        with open(os.path.join(tmp_runs_dir, "lruns_out"), "r") as out, open(
            os.path.join(tmp_runs_dir, "lruns_err"), "r"
        ) as err:
            stdout = out.read()
            stderr = err.read()
            return LRunsResult(stdout, stderr, returncode == 0)


@contextmanager
def runs_dir(request, runs, filenames, ignore_missing=False):
    "Copy run files into temporary directory, keeping the existing structure"
    with tempfile.TemporaryDirectory() as tmpdir:
        for run in runs:
            for filename in filenames:
                try:
                    with open_run_file(request, run, filename) as f:
                        rundir = os.path.join(tmpdir, run.root.root, run.number)
                        os.makedirs(rundir, exist_ok=True)
                        logger.debug(f"copying {run.filepath}/{filename} to {rundir}")
                        with open(os.path.join(rundir, filename), "wb") as fcopy:
                            shutil.copyfileobj(f, fcopy)
                except FileNotFoundError:
                    logger.debug(f"could not find {run.filepath}/{filename}")
                    if not ignore_missing:
                        raise
        yield tmpdir


def open_run_file(request, run: models.Run, filename: str):
    # Handle tutorial runs specially: all of their files are available within the app
    if FALSE and run.is_tutorial:
        return open(BASE_DIR / run.filepath / filename, "rb")
    else:
        data_product_uri = user_run_file_exists(request, run, filename)
        if data_product_uri is not None:
            return user_storage.open_file(request, data_product_uri=data_product_uri)
        else:
            raise FileNotFoundError(f"{filename} does not exist in run {run.id}")


def run_file_exists(request, run: models.Run, filename: str) -> bool:
    if FALSE and run.is_tutorial:
        return (BASE_DIR / run.filepath / filename).exists()
    else:
        return user_run_file_exists(request, run, filename) is not None


def user_run_file_exists(request, run, filename):
    """Return data product uri for run file if it exists, else None."""

    # check to see if the file is already in the run directory, for backwards compatibility
    if run.owner == request.user:
        print(f"DEBUG: Checking run directory for owner {run.owner}")
        try:
            data_product_uri = user_storage.user_file_exists(
                request, os.path.join(run.filepath, filename)
            )
            if data_product_uri is not None:
                 print(f"DEBUG: Found {filename} in {run.filepath} with URI: {data_product_uri}")
                 logger.debug(f"Found {filename} in {run.filepath}")
                 return data_product_uri
        except Exception as e:
            print(f"DEBUG: Failed to check run directory: {e}")
    else:
       print(f"DEBUG: Run owner {run.owner} != request user {request.user}, skipping run directory check")        

    # Find the most recent completed execution (Airavata experiment) if exists
    print(f"DEBUG: Looking for completed execution in {run.executions.count()} executions")
    experiment_model = None
    for execution in run.executions.order_by("-created"):
        status_name = execution.get_airavata_experiment_status(request)
        experiment_state = ExperimentState._NAMES_TO_VALUES[status_name]
        print(f"DEBUG: Execution {execution.airavata_experiment_id} status: {status_name}")
        if experiment_state == ExperimentState.COMPLETED:
            logger.debug(f"getExperiment({execution.airavata_experiment_id})")
            experiment_model = request.airavata_client.getExperiment(
                request.authz_token, execution.airavata_experiment_id
            )
            break
    if experiment_model is None:
        print(f"DEBUG: No completed execution found for run {run.id}")
        return None

    # Load the Modl_RunID file to find location of files
    # TODO: cache this information
    print(f"DEBUG: Looking for Modl_RunID in {len(experiment_model.experimentOutputs)} experiment outputs")
    modl_runid_output = None
    for output in experiment_model.experimentOutputs:
        if output.name == "Modl_RunID":
            modl_runid_output = output
            print(f"DEBUG: Found Modl_RunID output with URI: {output.value}")
            break
    if modl_runid_output is None:
        print(f"DEBUG: Modl_RunID output not found")
        raise Exception("Modl_RunID file is missing")
    if  not user_storage.exists(request, data_product_uri=modl_runid_output.value):
        print(f"DEBUG: Modl_RunID file does not exist at URI: {modl_runid_output.value}")
        raise Exception("Modl_RunID file is missing")

    print(f"DEBUG: Opening Modl_RunID file from URI: {modl_runid_output.value}")
    modl_runid_file = user_storage.open_file(
        request, data_product_uri=modl_runid_output.value
    )
    model_runid = modl_runid_file.read().decode()
    print(f"DEBUG: Modl_RunID file contents: '{model_runid.strip()}'")

    m = re.match(r"(\S+) (\S+)", model_runid)
    if m is None:
        print(f"DEBUG: Failed to parse Modl_RunID contents: '{model_runid}'")
        raise Exception(f"Invalid Modl_RunID file contents: {model_runid}")
    model, run_id = m.group(1, 2)
    print(f"DEBUG: Parsed model: '{model}', run_id: '{run_id}'")
    
    try:
        print(f"DEBUG: Using list_experiment_dir for experimen     t: '{experiment_model.experimentId}'")
        directories, files = user_storage.list_experiment_dir(     request, experiment_model.experimentId)
        print(f"DEBUG: Found {len(directories)} directories an     d {len(files)} files in experiment directory")
 
        available_files = [file['name'] for file in files]
        print(f"DEBUG: Available files: {available_files}") 
    # Check for the file in ARCHIVE/model/run_id/ directory
        for file in files:
            if file['name'] == filename:
                print(f"DEBUG: Found file '{filename}' with data product URI: '{file['data-product-uri']}'")
                return file['data-product-uri']
     
        print(f"DEBUG: File '{filename}' not found in experiment directory")
        return None 
         
    except Exception as e:
        print(f"DEBUG: Failed to list experiment directory: {e}")

        archive_path = os.path.join("ARCHIVE", model, run_id, filename)
        print(f"DEBUG: Fallback - checking archive path: '{archive_path}'")     
        
        try:
           data_product_uri = user_storage.user_file_exists(
               request,
               os.path.join("ARCHIVE", model, run_id, filename),
               experiment_id=experiment_model.experimentId,
           )
           return data_product_uri
        except Exception as fallback_e:
            print(f"DEBUG: Fallback also failed: {fallback_e}")
            return None


def get_run_output_data_product_uri(request, run: models.Run, data_type: str):
    # Find most recent execution
    if run.executions.exists():
        most_recent_execution: models.RemoteExecution = run.executions.order_by(
            "-created"
        )[0]
    else:
        return None

    experiment_model: ExperimentModel = request.airavata_client.getExperiment(
        request.authz_token, most_recent_execution.airavata_experiment_id
    )
    # Find the output by data type
    output = None
    for output in experiment_model.experimentOutputs:
        output_type_name = DataType._VALUES_TO_NAMES[output.type]
        if output_type_name == data_type:
            output = output
            break
    # If experiment is finished, see if there is an experimentOutput available
    if most_recent_execution.is_airavata_experiment_finished(request):
        if output is not None or user_storage.exists(
            request, data_product_uri=output.value
        ):
            return output.value
        else:
            return None
    # otherwise, see if there is an intermediate output available
    else:
        data_products = (
            experiment_util.intermediate_output.get_intermediate_output_data_products(
                request, experiment_model, output.name
            )
        )
        if len(data_products) == 1 and user_storage.exists(
            request, data_product=data_products[0]
        ):
            return data_products[0].productUri
        else:
            return None
