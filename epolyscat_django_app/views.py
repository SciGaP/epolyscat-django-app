"""REST endpoints for the ePolyScat app.

Runs and app experiments are stored server-side in Airavata via gRPC (see
run_store.py): a run is an ExperimentModel, a run-grouping "experiment" is an
Airavata Project. Only UI state (views, plot parameters, pointer rows) hits
the app's local database.
"""

import base64
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

from django.apps import apps
from django.conf import settings
from django.http import FileResponse
from django.shortcuts import render
from django.utils import timezone
from django_airavata.apps.auth.decorators import login_required
from rest_framework import (
    exceptions,
    pagination,
    permissions,
    response,
    viewsets,
)
from rest_framework.decorators import action, api_view
from rest_framework.response import Response

from epolyscat_django_app import (
    airavata_gateway,
    models,
    run_store,
    serializers,
    user_storage,
)

BASE_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
logger = logging.getLogger(__name__)


@login_required
def home(request):
    return render(
        request,
        "epolyscat_django_app/application.html",
        {"project_name": "ePolyScat-Django-App"},
    )


class IsOwner(permissions.BasePermission):
    def has_object_permission(self, request, view, obj):
        # Only the owner has write access
        return request.user.username == obj.owner


class IsOwnerOrReadOnly(permissions.BasePermission):
    def has_object_permission(self, request, view, obj):
        if request.method in permissions.SAFE_METHODS:
            return True
        # Only the owner has write access
        return request.user.username == obj.owner


class Pagination(pagination.PageNumberPagination):
    page_size = 10
    page_size_query_param = "page_size"


class ExperimentViewSet(viewsets.ViewSet):
    """Run-grouping "experiments", stored server-side as Airavata Projects.
    The local ProjectExtras registry marks which of the user's projects are
    app experiments (and tombstones deleted ones)."""

    serializer_class = serializers.ExperimentSerializer
    permission_classes = [permissions.IsAuthenticated]
    lookup_value_regex = "[^/]+"

    def _registered(self, request):
        """(project, extras) pairs for the user's app experiments, most
        recently touched first."""
        extras_by_id = {
            extras.project_id: extras
            for extras in models.ProjectExtras.objects.filter(deleted=False)
        }
        pairs = []
        for project in airavata_gateway.get_user_projects(request):
            project_id = airavata_gateway.project_id(project)
            if project_id in extras_by_id and project.owner == request.user.username:
                pairs.append((project, extras_by_id[project_id]))
        pairs.sort(key=lambda pair: pair[1].updated, reverse=True)
        return pairs

    def _get_pair(self, request, pk):
        extras = models.ProjectExtras.objects.filter(
            project_id=pk, deleted=False
        ).first()
        if extras is None:
            raise exceptions.NotFound(f"No experiment {pk}")
        try:
            project = airavata_gateway.get_project(request, pk)
        except Exception as e:
            raise exceptions.NotFound(f"No experiment {pk}") from e
        return project, extras

    def list(self, request):
        pairs = self._registered(request)
        paginator = Pagination()
        page = paginator.paginate_queryset(pairs, request, view=self)
        serializer = serializers.ExperimentSerializer(
            page, many=True, context={"request": request}
        )
        return paginator.get_paginated_response(serializer.data)

    def create(self, request):
        serializer = serializers.ExperimentSerializer(
            data=request.data, context={"request": request}
        )
        serializer.is_valid(raise_exception=True)
        serializer.save()
        return Response(serializer.data)

    def retrieve(self, request, pk=None):
        pair = self._get_pair(request, pk)
        serializer = serializers.ExperimentSerializer(
            pair, context={"request": request}
        )
        return Response(serializer.data)

    def update(self, request, pk=None, partial=False):
        pair = self._get_pair(request, pk)
        if pair[0].owner != request.user.username:
            raise exceptions.PermissionDenied(
                "You can only update an experiment that you own"
            )
        serializer = serializers.ExperimentSerializer(
            pair, data=request.data, partial=partial, context={"request": request}
        )
        serializer.is_valid(raise_exception=True)
        serializer.save()
        return Response(serializer.data)

    def partial_update(self, request, pk=None):
        return self.update(request, pk=pk, partial=True)

    def destroy(self, request, pk=None):
        project, extras = self._get_pair(request, pk)
        if project.owner != request.user.username:
            raise exceptions.PermissionDenied(
                "You can only delete an experiment that you own"
            )
        # Soft delete: the Airavata project (and the runs in it) stay
        # server-side; the experiment just stops being listed.
        extras.deleted = True
        extras.save()
        return Response(status=204)

    @action(detail=False, methods=["get"])
    def statistics(self, request):
        "Counts of undeleted experiments and their runs"
        pairs = self._registered(request)
        project_ids = {airavata_gateway.project_id(p) for p, _ in pairs}
        run_summaries = run_store.list_run_summaries(request)
        runs_count = sum(
            1 for summary in run_summaries if summary.project_id in project_ids
        )
        return Response(
            {"experiments_count": len(pairs), "runs_count": runs_count}
        )


class RunViewSet(viewsets.ViewSet):
    """Runs, stored server-side as Airavata experiments (see run_store)."""

    serializer_class = serializers.RunSerializer
    permission_classes = [permissions.IsAuthenticated]
    lookup_value_regex = "[^/]+"

    DATA_TYPE_TO_FILENAME = {
        "STDOUT": "stdout",
        "STDERR": "stderr",
    }

    def _serialize(self, request, run_or_runs, many=False):
        return serializers.RunSerializer(
            run_or_runs, many=many, context={"request": request, "view": self}
        ).data

    def get_object(self):
        request = self.request
        pk = self.kwargs.get("pk")
        try:
            run = run_store.get_run(request, pk)
        except Exception as e:
            raise exceptions.NotFound(f"No run {pk}") from e
        if run.extras.deleted:
            raise exceptions.NotFound(f"No run {pk}")
        return run

    def _require_owner(self, request, run):
        if run.owner != request.user.username:
            raise exceptions.PermissionDenied(
                "Only the owner can modify this run"
            )

    def list(self, request):
        view_id = request.query_params.get("viewId")
        project_id = request.query_params.get("experiment")

        if view_id:
            view = (
                models.View.filter_by_user(request).filter(id=view_id).first()
            )
            if view is None:
                raise exceptions.NotFound(f"No view {view_id}")
            if view.type == "unsubmitted":
                run_ids = [
                    summary.experiment_id
                    for summary in run_store.list_run_summaries(
                        request, project_id=project_id
                    )
                    if run_store.summary_status(summary) == "CREATED"
                ]
            else:
                run_ids = list(
                    view.runs.filter(deleted=False)
                    .order_by("-created")
                    .values_list("experiment_id", flat=True)
                )
        else:
            run_ids = [
                summary.experiment_id
                for summary in run_store.list_run_summaries(
                    request, project_id=project_id
                )
            ]

        paginator = Pagination()
        page_ids = paginator.paginate_queryset(run_ids, request, view=self)
        runs = []
        for run_id in page_ids:
            try:
                runs.append(run_store.get_run(request, run_id))
            except Exception:
                logger.warning("Could not load run %s", run_id, exc_info=True)
        return paginator.get_paginated_response(
            self._serialize(request, runs, many=True)
        )

    def retrieve(self, request, pk=None):
        run = self.get_object()
        return Response(self._serialize(request, run))

    def create(self, request):
        serializer = serializers.RunSerializer(
            data=request.data, context={"request": request}
        )
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        if not data.get("name"):
            raise exceptions.ValidationError({"name": ["This field is required."]})

        run = run_store.create_run(
            request,
            name=data["name"],
            description=data.get("description", ""),
            experiment_project_id=data.get("experiment"),
            group_resource_profile_id=data.get("group_resource_profile_id"),
            compute_resource_id=data.get("compute_resource_id"),
            queue_name=data.get("queue_name"),
            core_count=data.get("core_count"),
            node_count=data.get("node_count"),
            walltime_limit=data.get("walltime_limit"),
            total_physical_memory=data.get("total_physical_memory"),
            inputs_data=request.data.get("inputs_data", []),
        )

        for view_id in request.data.get("viewIds", []):
            if int(view_id) < 0:
                tutorial_view = models.View.tutorial_view()
                if tutorial_view is not None:
                    run.extras.views.add(tutorial_view)
            else:
                view = models.View.objects.filter(pk=int(view_id)).first()
                if view is not None:
                    run.extras.views.add(view)

        return Response(self._serialize(request, run))

    def update(self, request, pk=None, partial=False):
        run = self.get_object()
        self._require_owner(request, run)
        serializer = serializers.RunSerializer(
            data=request.data, partial=True, context={"request": request}
        )
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        run = run_store.update_run(
            request,
            run,
            name=data.get("name"),
            description=data.get("description"),
            group_resource_profile_id=data.get("group_resource_profile_id"),
            compute_resource_id=data.get("compute_resource_id"),
            queue_name=data.get("queue_name"),
            core_count=data.get("core_count"),
            node_count=data.get("node_count"),
            walltime_limit=data.get("walltime_limit"),
            total_physical_memory=data.get("total_physical_memory"),
            is_email_notification_on=data.get("is_email_notification_on"),
            inputs_data=request.data.get("inputs_data"),
        )
        return Response(self._serialize(request, run))

    def partial_update(self, request, pk=None):
        return self.update(request, pk=pk, partial=True)

    def destroy(self, request, pk=None):
        run = self.get_object()
        self._require_owner(request, run)
        delete_associated = request.GET.get("deleteAssociated", "true") != "false"
        run_store.delete_run(request, run, delete_associated=delete_associated)
        return Response(status=204)

    @action(detail=True, methods=["POST"])
    def clone(self, request, pk=None):
        run = self.get_object()
        data = self._serialize(request, run)
        for field in ("id", "created", "updated", "deleted", "executions"):
            data.pop(field, None)
        data["name"] = f"{data['name']} CLONE"
        data["status"] = run_store.STATUS_UNSUBMITTED
        data["job_status"] = run_store.STATUS_UNSUBMITTED
        return Response(data)

    @action(detail=True, methods=["post", "patch"])
    def submit(self, request, pk=None):
        run = self.get_object()
        self._require_owner(request, run)
        if request.data:
            serializer = serializers.RunSerializer(
                data=request.data, partial=True, context={"request": request}
            )
            serializer.is_valid(raise_exception=True)
            data = serializer.validated_data
            run = run_store.update_run(
                request,
                run,
                group_resource_profile_id=data.get("group_resource_profile_id"),
                compute_resource_id=data.get("compute_resource_id"),
                queue_name=data.get("queue_name"),
                core_count=data.get("core_count"),
                node_count=data.get("node_count"),
                walltime_limit=data.get("walltime_limit"),
                total_physical_memory=data.get("total_physical_memory"),
                inputs_data=request.data.get("inputs_data"),
            )

        # Only the compute resource is genuinely required to launch (and
        # create_run defaults it); unset numeric scheduling values coalesce
        # to 0 and the scheduler applies its own defaults — the UI leaves
        # e.g. total_physical_memory blank routinely.
        if not run.compute_resource_id:
            raise exceptions.ValidationError(
                {"compute_resource_id": ["must be provided for submission"]}
            )

        run_store.submit_run(request, run)
        return Response(self._serialize(request, run))

    @action(detail=True, methods=["post"])
    def resubmit(self, request, pk=None):
        run = self.get_object()
        self._require_owner(request, run)
        if request.data:
            serializer = serializers.RunSerializer(
                data=request.data, partial=True, context={"request": request}
            )
            serializer.is_valid(raise_exception=True)
            data = serializer.validated_data
            run = run_store.update_run(
                request,
                run,
                group_resource_profile_id=data.get("group_resource_profile_id"),
                compute_resource_id=data.get("compute_resource_id"),
                queue_name=data.get("queue_name"),
                core_count=data.get("core_count"),
                node_count=data.get("node_count"),
                walltime_limit=data.get("walltime_limit"),
                total_physical_memory=data.get("total_physical_memory"),
            )
        run_store.resubmit_run(request, run)
        return Response(self._serialize(request, run))

    @action(detail=True, methods=["put"])
    def cancel(self, request, pk=None):
        run = self.get_object()
        self._require_owner(request, run)
        run_store.cancel_run(request, run)
        return Response(self._serialize(request, run))

    @action(detail=True, methods=["get"])
    def status(self, request, pk=None):
        run = self.get_object()
        return Response(run.status())

    @action(detail=True, methods=["PATCH"])
    def change_notification_settings(self, request, pk=None):
        run = self.get_object()
        self._require_owner(request, run)
        run = run_store.update_run(
            request,
            run,
            is_email_notification_on=request.data.get(
                "is_email_notification_on", False
            ),
        )
        return Response(self._serialize(request, run))

    @action(detail=False, methods=["GET"])
    def tutorial_runs(self, request):
        tutorial_view = models.View.tutorial_view()
        if tutorial_view is None:
            return Response([])
        runs = []
        for extras in tutorial_view.runs.filter(deleted=False):
            try:
                run = run_store.get_run(request, extras.experiment_id)
                run.extras = extras
                runs.append(run)
            except Exception:
                logger.warning(
                    "Could not load tutorial run %s", extras.experiment_id
                )
        return Response(self._serialize(request, runs, many=True))

    @action(methods=["get"], detail=True)
    def viewables(self, request, pk=None):
        run = self.get_object()
        viewables = []

        epolyscat_settings = apps.get_app_config(
            "epolyscat_django_app"
        ).APPLICATION_SETTINGS["EPOLYSCAT_DJANGO_APP"]
        for filename, _description in epolyscat_settings["FILE_VIEWABLE"].items():
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
            except Exception:
                logger.exception(f"Failed to check if {data_type} is available")

        return Response(viewables)

    @action(methods=["get"], detail=True, url_path="input-files")
    def input_files(self, request, pk=None):
        run = self.get_object()
        input_files_list = []

        epolyscat_settings = apps.get_app_config(
            "epolyscat_django_app"
        ).APPLICATION_SETTINGS["EPOLYSCAT_DJANGO_APP"]
        for filename, _description in epolyscat_settings["FILE_INPUT"].items():
            if run_file_exists(request, run, filename):
                url = self.reverse_action(url_name="show-viewable", args=[pk, filename])
                input_files_list.append(dict(filename=filename, url=url))
        return Response(input_files_list)

    @action(methods=["get"], detail=True, url_path=r"viewables/(?P<filename>[\w]+)")
    def show_viewable(self, request, pk=None, filename: str | None = None):
        run = self.get_object()
        # The url_path regex always captures a filename; the default is only
        # here to match the action signature.
        if filename is None:
            raise exceptions.NotFound("No filename given")
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
            return FileResponse(f, content_type="text/plain; charset=utf-8")
        except FileNotFoundError as e:
            raise exceptions.NotFound from e

    @action(detail=True, methods=["GET"])
    def get_output_files(self, request, pk=None):
        run = self.get_object()
        output_files = []

        latest = run.latest_execution_id
        if latest is not None:
            try:
                output_files = user_storage.list_experiment_dir(
                    request, latest, path="ARCHIVE"
                )[1]
            except Exception:
                pass

        return Response(output_files)


class PlotParametersViewSet(viewsets.ModelViewSet):
    serializer_class = serializers.PlotParametersSerializer
    permission_classes = [permissions.IsAuthenticated, IsOwner]

    def get_queryset(self):
        request = self.request
        return models.PlotParameters.filter_by_user(request).order_by("-last_use")


class ViewsViewSet(viewsets.ModelViewSet):
    serializer_class = serializers.ViewSerializer
    permission_classes = [permissions.IsAuthenticated, IsOwnerOrReadOnly]
    pagination_class = Pagination

    def get_queryset(self):
        return models.View.objects.filter(
            owner=self.request.user.username, deleted=False
        ).order_by("-order", "-updated")

    def _add_runs(self, request, view, run_ids):
        for run_id in run_ids:
            extras = models.RunExtras.for_experiment(run_id)
            extras.views.add(view)

    def create(self, request, *args, **kwargs):
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        view = serializer.save(owner=request.user.username)
        self._add_runs(request, view, request.data.get("runIds", []))
        return Response(self.get_serializer(view).data)

    def update(self, request, *args, **kwargs):
        kwargs.pop("partial", False)
        view = self.get_object()
        serializer = self.get_serializer(view, data=request.data, partial=True)
        serializer.is_valid(raise_exception=True)
        serializer.save()

        if request.data.get("overide"):
            view.runs.clear()
        self._add_runs(request, view, request.data.get("runIds", []))
        # Membership changed: rebuild the cached server runs on next access.
        if hasattr(view, "_server_runs"):
            del view._server_runs
        return Response(self.get_serializer(view).data)

    def perform_destroy(self, instance):
        if instance.type != "user-defined":
            raise exceptions.ValidationError("Can only delete 'user-defined' views.")
        instance.deleted = True
        instance.save()

    @action(methods=["GET"], detail=False)
    def tutorials(self, request):
        tutorials_view = models.View.tutorial_view()
        if tutorials_view is None:
            tutorials_view = models.View.objects.create(
                name="Tutorials", type="tutorial", owner=None
            )
        serializer = self.get_serializer(tutorials_view)
        return Response(serializer.data)


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

    plotfiles = serializer.validated_data["plotfiles"]

    def map_plotfile_to_obj(plotfile):
        file = user_storage.open_file(
            request, data_product_uri=plotfile["data_product_uri"]
        )
        sep = "" if plotfile["prefix"] == "" else "__"
        file_name = file.name.split("/")[-1]
        return {"file": file, "name": f"{plotfile['prefix']}{sep}{file_name}"}

    plotfile_objs = list(map(map_plotfile_to_obj, plotfiles))
    with get_tmp_dir(request, plotfile_objs, False) as tmp_dir:
        plot_parameters = serializer.validated_data.get("plot_parameters", None)
        if plot_parameters is not None:
            obj, created = models.PlotParameters.objects.get_or_create(
                xaxis=plot_parameters["xaxis"],
                yaxes=plot_parameters["yaxes"],
                flags=plot_parameters["flags"],
                owner=request.user.username,
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
                    cwd=tmp_dir,
                )
                process.wait()
                process_output = process.stdout.read().decode()
                logger.debug(f"plot.py -showColumns returncode={process.returncode}")
                logger.debug(f"plot.py -showColumns output={process_output}")
                data["user_guidance"] = process_output

            return response.Response(data)
        except Exception:
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


def open_run_file(request, run: run_store.ServerRun, filename: str):
    data_product_uri = user_run_file_exists(request, run, filename)
    if data_product_uri is not None:
        return user_storage.open_file(request, data_product_uri=data_product_uri)
    else:
        raise FileNotFoundError(f"{filename} does not exist in run {run.id}")


def run_file_exists(request, run: run_store.ServerRun, filename: str) -> bool:
    return user_run_file_exists(request, run, filename) is not None


def user_run_file_exists(request, run: run_store.ServerRun, filename):
    """Return data product uri for run file if it exists, else None."""

    # Check the run's own storage directory first.
    if run.owner == request.user.username:
        try:
            data_product_uri = user_storage.user_file_exists(
                request, os.path.join(run.directory, filename)
            )
            if data_product_uri is not None:
                logger.debug(f"Found {filename} in {run.directory}")
                return data_product_uri
        except Exception:
            logger.debug("Failed to check run directory", exc_info=True)

    # Find the most recent completed execution (Airavata experiment) if any.
    experiment_model = None
    for execution_id in reversed(run.execution_ids):
        current = airavata_gateway.get_experiment_status(request, execution_id)
        if current.state == airavata_gateway.experiment_completed_state():
            experiment_model = airavata_gateway.get_experiment(request, execution_id)
            break
    if experiment_model is None:
        return None

    # Load the Modl_RunID file to find the location of output files.
    modl_runid_output = None
    for output in airavata_gateway.experiment_outputs(experiment_model):
        if output.name == "Modl_RunID":
            modl_runid_output = output
            break
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

    experiment_id = airavata_gateway.experiment_id(experiment_model)
    try:
        _directories, files = user_storage.list_experiment_dir(
            request, experiment_id
        )
        for file in files:
            if file["name"] == filename:
                return file.get("data-product-uri") or file.get("data_product_uri")
        return None
    except Exception:
        logger.debug("Failed to list experiment directory", exc_info=True)
        try:
            return user_storage.user_file_exists(
                request,
                os.path.join("ARCHIVE", model, run_id, filename),
                experiment_id=experiment_id,
            )
        except Exception:
            logger.debug("Archive fallback failed", exc_info=True)
            return None


def get_run_output_data_product_uri(request, run: run_store.ServerRun, data_type: str):
    latest = run.latest_execution_id
    if latest is None:
        return None

    experiment_model = airavata_gateway.get_experiment(request, latest)
    # Find the output by data type
    output = None
    for candidate in airavata_gateway.experiment_outputs(experiment_model):
        if airavata_gateway.data_type_name(candidate.type) == data_type:
            output = candidate
            break
    if output is None:
        return None
    # If experiment is finished, see if there is an experimentOutput available
    if run.is_execution_finished(latest):
        if output.value and user_storage.exists(
            request, data_product_uri=output.value
        ):
            return output.value
        return None
    # otherwise, see if there is an intermediate output available
    data_products = airavata_gateway.get_intermediate_output_data_products(
        request, experiment_model, output.name
    )
    if len(data_products) == 1 and user_storage.exists(
        request,
        data_product_uri=airavata_gateway.data_product_uri(data_products[0]),
    ):
        return airavata_gateway.data_product_uri(data_products[0])
    return None


@api_view(["GET"])
def api_settings(request):
    app_module_id = getattr(settings, "EPOLYSCAT", {}).get(
        "EPOLYSCAT_APPLICATION_ID",
        "ePolyScat_940ab1c9-4ceb-431c-8595-c6246a195442",
    )
    return response.Response({"EPOLYSCAT": {"EPOLYSCAT_APPLICATION_ID": app_module_id}})
