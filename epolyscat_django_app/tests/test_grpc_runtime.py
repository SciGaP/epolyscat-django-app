from types import SimpleNamespace
from unittest import mock

from django.contrib.auth import get_user_model
from django.test import RequestFactory, TestCase, override_settings

from epolyscat_django_app import models, serializers, views


class _CopyTarget:
    __slots__ = ("value",)

    def __init__(self):
        self.value = None

    def CopyFrom(self, value):
        self.value = value


class _GrpcExperiment:
    __slots__ = (
        "experiment_name",
        "execution_id",
        "gateway_id",
        "user_name",
        "project_id",
        "experiment_inputs",
        "experiment_outputs",
        "user_configuration_data",
    )

    def __init__(
        self,
        experiment_name="",
        execution_id="",
        gateway_id="",
        user_name="",
    ):
        self.experiment_name = experiment_name
        self.execution_id = execution_id
        self.gateway_id = gateway_id
        self.user_name = user_name
        self.project_id = ""
        self.experiment_inputs = []
        self.experiment_outputs = []
        self.user_configuration_data = _CopyTarget()


class _GrpcInput:
    __slots__ = (
        "name",
        "type",
        "value",
        "is_required",
        "required_to_added_to_command_line",
        "input_order",
    )

    def __init__(self, name, input_type, input_order):
        self.name = name
        self.type = input_type
        self.value = ""
        self.is_required = True
        self.required_to_added_to_command_line = False
        self.input_order = input_order


class _GrpcUserConfiguration:
    __slots__ = (
        "group_resource_profile_id",
        "share_experiment_publicly",
        "computational_resource_scheduling",
    )

    def __init__(
        self,
        group_resource_profile_id,
        share_experiment_publicly,
        computational_resource_scheduling,
    ):
        self.group_resource_profile_id = group_resource_profile_id
        self.share_experiment_publicly = share_experiment_publicly
        self.computational_resource_scheduling = computational_resource_scheduling


class _GrpcScheduling:
    __slots__ = (
        "resource_host_id",
        "total_cpu_count",
        "node_count",
        "wall_time_limit",
        "queue_name",
        "total_physical_memory",
    )

    def __init__(self, **kwargs):
        for name in self.__slots__:
            setattr(self, name, kwargs.get(name))


class _GrpcProject:
    __slots__ = ("owner", "gateway_id", "name", "project_id")

    def __init__(
        self,
        owner="",
        gateway_id="",
        name="",
        project_id="",
    ):
        self.owner = owner
        self.gateway_id = gateway_id
        self.name = name
        self.project_id = project_id


class _GrpcJob:
    __slots__ = ("job_id",)

    def __init__(self, job_id):
        self.job_id = job_id


class _GrpcJobStatus:
    __slots__ = ("job_state", "time_of_state_change")

    def __init__(self, job_state, time_of_state_change):
        self.job_state = job_state
        self.time_of_state_change = time_of_state_change


class _GrpcJobStatusesResponse:
    __slots__ = ("job_statuses",)

    def __init__(self, job_statuses):
        self.job_statuses = job_statuses


class GrpcRuntimeTests(TestCase):
    def _create_run(self, user):
        experiment = models.Experiment.objects.create(
            name="experiment",
            description="",
            owner=user,
            airavata_project_id="project-id",
        )
        root = models.RunsRoot.objects.create(root="root", owner=user)
        return models.Run.objects.create(
            name="run",
            owner=user,
            root=root,
            number="0001",
            experiment=experiment,
            group_resource_profile_id="group-id",
            compute_resource_id="compute-id",
            queue_name="normal",
            core_count=8,
            node_count=2,
            walltime_limit=120,
            total_physical_memory=4096,
        )

    @override_settings(GATEWAY_ID="test-gateway")
    @mock.patch.object(views, "ExperimentModel", _GrpcExperiment)
    @mock.patch.object(views, "UserConfigurationDataModel", _GrpcUserConfiguration)
    @mock.patch.object(
        views,
        "ComputationalResourceSchedulingModel",
        _GrpcScheduling,
    )
    def test_remote_execution_uses_grpc_facade_and_proto_fields(self):
        user = get_user_model().objects.create_user(username="grpc-user")
        request = RequestFactory().post(
            "/epolyscat_django_app/api/runs/1/submit/"
        )
        request.user = user
        inputs = [_GrpcInput("Input_File", views.DataType.URI, 0)]
        research = SimpleNamespace(
            get_application_interface=mock.Mock(
                return_value=SimpleNamespace(
                    application_inputs=inputs,
                    application_outputs=[],
                )
            ),
            create_experiment=mock.Mock(return_value="experiment-id"),
            launch_experiment=mock.Mock(),
        )
        compute = SimpleNamespace(
            get_compute_resource=mock.Mock(
                return_value=SimpleNamespace(host_name="cluster.example.org")
            )
        )
        request.airavata = SimpleNamespace(research=research, compute=compute)
        run = self._create_run(user)

        execution = views.RunViewSet()._create_remote_execution(
            request,
            run,
            {"Input_File": "airavata-dp://input"},
            "app-interface-id",
            {"Input_File": "airavata-dp://input"},
            is_tutorial=False,
        )

        created_experiment = research.create_experiment.call_args.args[1]
        self.assertEqual(created_experiment.execution_id, "app-interface-id")
        self.assertEqual(created_experiment.gateway_id, "test-gateway")
        self.assertEqual(created_experiment.project_id, "project-id")
        self.assertEqual(
            created_experiment.experiment_inputs[0].value,
            "airavata-dp://input",
        )
        scheduling = (
            created_experiment.user_configuration_data.value
            .computational_resource_scheduling
        )
        self.assertEqual(scheduling.total_physical_memory, 4096)
        self.assertEqual(execution.resource_name, "cluster.example.org")
        research.get_application_interface.assert_called_once_with(
            "app-interface-id"
        )
        research.create_experiment.assert_called_once_with(
            "test-gateway",
            created_experiment,
        )
        research.launch_experiment.assert_called_once_with(
            "experiment-id",
            "test-gateway",
        )
        compute.get_compute_resource.assert_called_once_with("compute-id")

    @override_settings(GATEWAY_ID="test-gateway")
    @mock.patch.object(models, "AiravataProject", _GrpcProject)
    def test_experiment_project_creation_uses_grpc_facade_and_proto_fields(self):
        user = get_user_model().objects.create_user(username="project-user")
        research = SimpleNamespace(
            create_project=mock.Mock(return_value="project-id")
        )
        request = RequestFactory().post("/epolyscat_django_app/api/runs/")
        request.user = user
        request.airavata = SimpleNamespace(research=research)
        experiment = models.Experiment(name="science", owner=user)

        experiment.create_airavata_project(request)

        created_project = research.create_project.call_args.args[1]
        self.assertEqual(created_project.gateway_id, "test-gateway")
        self.assertEqual(created_project.owner, "project-user")
        self.assertEqual(experiment.airavata_project_id, "project-id")
        research.create_project.assert_called_once_with(
            "test-gateway",
            created_project,
        )

    @override_settings(GATEWAY_ID="test-gateway")
    @mock.patch.object(serializers, "Project", _GrpcProject)
    def test_run_creation_selects_grpc_project_by_proto_field(self):
        user = get_user_model().objects.create_user(username="run-project-user")
        research = SimpleNamespace(
            get_user_projects=mock.Mock(
                return_value=[
                    _GrpcProject(project_id="Default_Project_existing")
                ]
            ),
            create_project=mock.Mock(),
        )
        request = RequestFactory().post("/epolyscat_django_app/api/runs/")
        request.user = user
        request.airavata = SimpleNamespace(research=research)

        run = serializers.RunSerializer(
            context={"request": request}
        ).create({"name": "grpc run", "owner": user})

        self.assertEqual(run.airavata_project_id, "Default_Project_existing")
        research.get_user_projects.assert_called_once_with(
            "test-gateway",
            "run-project-user",
            -1,
            0,
        )
        research.create_project.assert_not_called()

    @override_settings(GATEWAY_ID="test-gateway")
    def test_remote_execution_status_job_and_cancel_use_grpc_facade(self):
        user = get_user_model().objects.create_user(username="status-user")
        run = self._create_run(user)
        execution = models.RemoteExecution.objects.create(
            run=run,
            airavata_experiment_id="experiment-id",
        )
        research = SimpleNamespace(
            get_experiment_status=mock.Mock(
                return_value=SimpleNamespace(
                    state=models.ExperimentState.COMPLETED
                )
            ),
            get_job_details=mock.Mock(
                return_value=[_GrpcJob("job-id")]
            ),
            terminate_experiment=mock.Mock(),
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/")
        request.user = user
        request.airavata = SimpleNamespace(research=research)

        self.assertEqual(
            execution.get_airavata_experiment_status(request),
            "COMPLETED",
        )
        self.assertEqual(execution.get_job_id(request), "job-id")
        execution.cancel(request)

        research.get_experiment_status.assert_called_once_with("experiment-id")
        research.get_job_details.assert_called_once_with("experiment-id")
        research.terminate_experiment.assert_called_once_with(
            "experiment-id",
            "test-gateway",
        )

    def test_run_serializer_reads_job_status_map_from_grpc_response(self):
        user = get_user_model().objects.create_user(username="job-status-user")
        run = self._create_run(user)
        models.RemoteExecution.objects.create(
            run=run,
            airavata_experiment_id="experiment-id",
            airavata_experiment_status="COMPLETED",
        )
        research = SimpleNamespace(
            get_job_statuses=mock.Mock(
                return_value=_GrpcJobStatusesResponse(
                    {
                        "job-id": _GrpcJobStatus(
                            job_state=3,
                            time_of_state_change=20,
                        )
                    }
                )
            )
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/")
        request.user = user
        request.airavata = SimpleNamespace(research=research)

        status = serializers.RunSerializer(
            context={"request": request}
        ).get_job_status(run)

        self.assertEqual(status, "COMPLETED")
        research.get_job_statuses.assert_called_once_with("experiment-id")

    @override_settings(GATEWAY_ID="test-gateway")
    def test_runtime_catalog_queries_use_grpc_research_facade(self):
        research = SimpleNamespace(
            get_all_app_modules=mock.Mock(return_value=[]),
            get_all_application_interfaces=mock.Mock(
                return_value=[
                    SimpleNamespace(
                        application_interface_id="interface-id",
                        application_modules=["module-id"],
                    )
                ]
            ),
            get_all_application_deployments=mock.Mock(
                return_value=[
                    SimpleNamespace(
                        app_module_id="module-id",
                        compute_host_id="compute-id",
                        executable_path="/apps/module/run.sh",
                    )
                ]
            ),
        )
        request = RequestFactory().get(
            "/epolyscat_django_app/api/runs/runtime_audit/"
        )
        request.airavata = SimpleNamespace(research=research)
        viewset = views.RunViewSet()

        self.assertEqual(
            viewset._get_app_interface_id_for_module(request, "module-id"),
            "interface-id",
        )
        run = SimpleNamespace(compute_resource_id="compute-id")
        with mock.patch.object(
            viewset,
            "_get_run_application_module_id",
            return_value="module-id",
        ):
            self.assertEqual(
                viewset._get_run_deployment_executable_path(request, run),
                "/apps/module/run.sh",
            )
        response = viewset.runtime_audit(request)

        self.assertIn("ready", response.data)
        research.get_all_app_modules.assert_called_once_with("test-gateway")
        self.assertGreaterEqual(
            research.get_all_application_interfaces.call_count,
            2,
        )
        self.assertGreaterEqual(
            research.get_all_application_deployments.call_count,
            2,
        )

    def test_output_lookup_uses_grpc_experiment_fields(self):
        user = get_user_model().objects.create_user(username="output-user")
        run = self._create_run(user)
        models.RemoteExecution.objects.create(
            run=run,
            airavata_experiment_id="experiment-id",
            airavata_experiment_status="COMPLETED",
        )
        output = SimpleNamespace(
            name="Output",
            type=views.DataType.URI,
            value="airavata-dp://output",
        )
        research = SimpleNamespace(
            get_experiment=mock.Mock(
                return_value=SimpleNamespace(experiment_outputs=[output])
            )
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/")
        request.user = user
        request.airavata = SimpleNamespace(research=research)

        with mock.patch.object(
            views.user_storage,
            "exists",
            return_value=True,
        ):
            product_uri = views.get_run_output_data_product_uri(
                request,
                run,
                "URI",
            )

        self.assertEqual(product_uri, "airavata-dp://output")
        research.get_experiment.assert_called_once_with("experiment-id")
