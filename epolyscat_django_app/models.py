import logging
import re
from typing import Union

from airavata.model.experiment.ttypes import ExperimentModel
from airavata.model.status.ttypes import ExperimentState
from airavata.model.workspace.ttypes import Project
from airavata_django_portal_sdk import experiment_util, user_storage
from django.conf import settings
from django.db import models
from django.db.models import Q

logger = logging.getLogger(__name__)


class Experiment(models.Model):
    name = models.CharField(max_length=255)
    description = models.CharField(max_length=4000)
    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.CASCADE,
        null=True,
        related_name="epolyscat_experiments"
    )
    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)
    deleted = models.BooleanField(default=False)
    airavata_project_id = models.CharField(max_length=255, unique=True, null=True)
    root = models.OneToOneField(
        "RunsRoot", on_delete=models.CASCADE, related_name="experiment"
    )

    class Meta:
        unique_together = ["name", "owner"]

    def create_airavata_project(self, request):
        airavata_project = Project(
            owner=request.user.username,
            gatewayId=settings.GATEWAY_ID,
            name="Runs for experiment: " + self.name,
        )
        airavata_project_id = request.airavata_client.createProject(
            request.authz_token, settings.GATEWAY_ID, airavata_project
        )
        self.airavata_project_id = airavata_project_id

    def __str__(self):
        return self.name


class RunsRoot(models.Model):
    root = models.CharField(max_length=100)
    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.CASCADE,
        null=True,
        related_name="epolyscat_runsroots"
    )
    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)

    class Meta:
        unique_together = ["root", "owner"]

    def get_next_run_number(self):
        max_number = max(
            map(lambda r: int(r["number"]), self.runs.values("number")), default=0
        )
        return f"{max_number+1:04}"

    def __str__(self):
        return self.root


class View(models.Model):
    name = models.CharField(max_length=255)
    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.CASCADE,
        null=True,
        related_name="epolyscat_views"
    )
    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)
    VIEW_TYPE_CHOICES = (
        ("unsubmitted", "Unsubmitted"),
        ("default", "Selected"),
        ("user-defined", "User Defined"),
        ("tutorial", "Tutorial"),
    )
    type = models.CharField(
        max_length=20, default="user-defined", choices=VIEW_TYPE_CHOICES
    )
    deleted = models.BooleanField(default=False)
    order = models.IntegerField(default=0)

    class Meta:
        unique_together = ["name", "owner"]

    def __str__(self):
        return self.name

    @staticmethod
    def filter_by_user(request):
        return View.objects.filter(Q(owner=request.user) | Q(owner=None))

    def populate_unsubmitted_runs(self, request):
        executions = RemoteExecution.objects.filter(run=models.OuterRef("pk"))
        self.runs.set(Run.filter_by_user(request).filter(~models.Exists(executions)))

    @staticmethod
    def create_default_views(request):
        owner = request.user
        if not View.objects.filter(owner=owner, type="unsubmitted").exists():
            View.objects.create(
                type="unsubmitted", name="Unsubmitted", owner=owner, order=20
            )
        if not View.objects.filter(owner=owner, type="default").exists():
            View.objects.create(type="default", name="Selected", owner=owner, order=10)


class Run(models.Model):
    number = models.CharField(max_length=64)
    root = models.ForeignKey(RunsRoot, on_delete=models.CASCADE, related_name="runs")
    experiment = models.ForeignKey(
        Experiment, on_delete=models.CASCADE, related_name="runs"
    )
    views = models.ManyToManyField(View, related_name="runs")
    filepath = models.CharField(max_length=255)
    inpc_data_product_uri = models.CharField(max_length=64, null=True, default=None)
    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)
    deleted = models.BooleanField(default=False)
    group_resource_profile_id = models.CharField(max_length=255, null=True)
    compute_resource_id = models.CharField(max_length=255, null=True)
    queue_name = models.CharField(max_length=64, null=True)
    core_count = models.IntegerField(null=True)
    node_count = models.IntegerField(null=True)
    walltime_limit = models.IntegerField(null=True)
    total_physical_memory = models.IntegerField(null=True)
    input_table = models.TextField(null=True)

    class Meta:
        unique_together = ["number", "root"]
        ordering = ["root__root", "number"]

    @property
    def owner(self):
        return self.experiment.owner

    @property
    def name(self):
        return f"{self.root.root}/{self.number}"

    @property
    def is_tutorial(self):
        return Run.check_is_tutorial(self.id)

    @staticmethod
    def check_is_tutorial(run_id):
        try:
            tutorials_view = View.objects.get(type="tutorial", owner=None)
            return tutorials_view.runs.filter(pk=run_id).exists()
        except View.DoesNotExist:
            return False

    @staticmethod
    def filter_by_user(request):
        return Run.objects.filter(
            Q(experiment__owner=request.user)
            | Q(experiment__owner=None)
        )

    def get_most_recent_job_id(self, request):
        job_id = None
        for execution in self.executions.order_by("-created"):
            job_id = execution.get_job_id(request)
            if job_id is not None:
                break
        return job_id

    def are_all_executions_finished(self, request):
        for execution in self.executions.all():
            if not execution.is_airavata_experiment_finished(request):
                return False
        return True

    @property
    def latest_execution(self):
        try:
            return self.executions.latest("created")
        except RemoteExecution.DoesNotExist:
            return None

    def is_cancelable(self, request):
        latest_execution: RemoteExecution = self.latest_execution
        return latest_execution is not None and latest_execution.is_cancelable(request)

    def __str__(self):
        return self.number


class RemoteExecution(models.Model):
    run = models.ForeignKey(Run, on_delete=models.CASCADE, related_name="executions")
    airavata_experiment_id = models.CharField(max_length=255, unique=True, null=True)
    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)
    airavata_experiment_status = models.CharField(
        max_length=255,
        default=ExperimentState._VALUES_TO_NAMES[ExperimentState.CREATED],
    )
    resource_name = models.CharField(max_length=255, blank=True, default="")
    job_id = models.CharField(max_length=255, null=True)

    def get_job_id(self, request):
        if self.job_id is None and self.airavata_experiment_id is not None:
            jobs = request.airavata_client.getJobDetails(
                request.authz_token, self.airavata_experiment_id
            )
            if len(jobs) > 0:
                self.job_id = jobs[0].jobId
                self.save()
        return self.job_id

    def get_airavata_experiment_status(self, request):
        terminal_states = self.get_airavata_experiment_terminal_states()
        old_state = ExperimentState._NAMES_TO_VALUES[self.airavata_experiment_status]
        if old_state in terminal_states:
            return self.airavata_experiment_status
        else:
            logger.debug(f"getExperimentStatus({self.airavata_experiment_id})")
            current_status = request.airavata_client.getExperimentStatus(
                request.authz_token, self.airavata_experiment_id
            )
            self.airavata_experiment_status = ExperimentState._VALUES_TO_NAMES[
                current_status.state
            ]
            self.save()
            return self.airavata_experiment_status

    def is_airavata_experiment_finished(self, request):
        status = self.get_airavata_experiment_status(request)
        state = ExperimentState._NAMES_TO_VALUES[status]
        return state in self.get_airavata_experiment_terminal_states()

    def get_application_specific_status(self, request) -> Union[str, None]:
        output_name = "tRecX_Standard-Out"
        experiment_model: ExperimentModel = request.airavata_client.getExperiment(
            request.authz_token, self.airavata_experiment_id
        )
        stdout_dp = None
        if experiment_model.experimentStatus[-1].state == ExperimentState.COMPLETED:
            for output in experiment_model.experimentOutputs:
                if (
                    output.name == output_name
                    and output.value
                    and output.value.startswith("airavata-dp://")
                ):
                    stdout_dp = request.airavata_client.getDataProduct(
                        request.authz_token, output.value
                    )
        else:
            can_fetch = (
                experiment_util.intermediate_output.can_fetch_intermediate_output(
                    request, experiment_model, output_name
                )
            )
            if can_fetch:
                experiment_util.intermediate_output.fetch_intermediate_output(
                    request, self.airavata_experiment_id, output_name
                )
            data_products = experiment_util.intermediate_output.get_intermediate_output_data_products(
                request, experiment_model, output_name
            )
            if len(data_products) == 1:
                stdout_dp = data_products[0]

        if stdout_dp is not None:
            with user_storage.open_file(request, data_product=stdout_dp) as f:
                stdout = list(map(bytes.decode, f.readlines()))

                ident = re.compile(r"^ ===  done")
                items = list(filter(ident.match, stdout))
                if items or stdout[-1].find(" *** further output on") != -1:
                    return "* - Completed"

                ident = re.compile(r"^ ===  TIME PROPAGATION")
                items = list(filter(ident.match, stdout))
                if items:
                    return "Time propagating"

                ident = re.compile(r"^ ===  OUTPUT")
                items = list(filter(ident.match, stdout))
                if items:
                    return "Started"

        return None

    def get_airavata_experiment_terminal_states(self):
        return [
            ExperimentState.CANCELED,
            ExperimentState.COMPLETED,
            ExperimentState.FAILED,
        ]

    def get_airavata_experiment_cancelable_states(self):
        return [
            ExperimentState.VALIDATED,
            ExperimentState.SCHEDULED,
            ExperimentState.LAUNCHED,
            ExperimentState.EXECUTING,
        ]

    @property
    def resource_name_short(self):
        return self.resource_name.split(".")[0]

    def cancel(self, request):
        request.airavata_client.terminateExperiment(
            request.authz_token, self.airavata_experiment_id, settings.GATEWAY_ID
        )

    def is_cancelable(self, request) -> bool:
        status = self.get_airavata_experiment_status(request)
        state = ExperimentState._NAMES_TO_VALUES[status]
        return state in self.get_airavata_experiment_cancelable_states()


class PlotParameters(models.Model):
    xaxis = models.CharField(max_length=2, default="")
    yaxes = models.CharField(max_length=20, default="")
    flags = models.CharField(max_length=200, default="")
    last_use = models.DateTimeField(auto_now=True)
    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.CASCADE,
        related_name="epolyscat_plotparameters"
    )
    created = models.DateTimeField(auto_now_add=True)

    @staticmethod
    def filter_by_user(request):
        return PlotParameters.objects.filter(owner=request.user)

    def __str__(self) -> str:
        return f"x={self.xaxis} y={self.yaxes} {self.flags}"
