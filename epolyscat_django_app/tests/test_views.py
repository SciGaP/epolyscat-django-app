from io import StringIO
from types import SimpleNamespace
from unittest import mock

from django.contrib.auth import get_user_model
from django.test import RequestFactory, TestCase, override_settings
from rest_framework import exceptions

from epolyscat_django_app import models, serializers, views


class RunViewSetBackendTests(TestCase):
    def create_run(self, user):
        experiment = models.Experiment.objects.create(
            name="experiment",
            description="",
            owner=user,
        )
        root = models.RunsRoot.objects.create(root="root", owner=user)
        return models.Run.objects.create(
            name="run",
            owner=user,
            root=root,
            number="0001",
            experiment=experiment,
        )

    def test_run_viewset_exposes_workflow_continuation_action(self):
        action_names = {
            action.__name__ for action in views.RunViewSet.get_extra_actions()
        }

        self.assertIn("workflow_continuation", action_names)

    def test_workflow_continuation_get_returns_backend_eligibility(self):
        user = get_user_model().objects.create_user(username="continue-get")
        source = self.create_run(user)
        source.module_application = "Gaussian16"
        source.save(update_fields=["module_application"])
        models.RemoteExecution.objects.create(
            run=source,
            airavata_experiment_id="completed-gaussian",
            airavata_experiment_status="COMPLETED",
        )
        request = RequestFactory().get(
            f"/epolyscat_django_app/api/runs/{source.id}/workflow_continuation/"
        )
        request.user = user
        viewset = views.RunViewSet()
        viewset.request = request
        viewset.get_object = mock.Mock(return_value=source)

        response = viewset.workflow_continuation(request, pk=source.id)

        self.assertIn("eligible", response.data)
        self.assertTrue(response.data["eligible"])
        self.assertEqual(response.data["source_stage"], "Data_Gen")
        self.assertEqual(response.data["next_stage"], "ePolyScat_Run")

    def test_workflow_continuation_post_creates_parent_after_imported_source(self):
        user = get_user_model().objects.create_user(username="continue-post")
        source = self.create_run(user)
        source.module_application = "ePolyScat"
        source.group_resource_profile_id = "group"
        source.compute_resource_id = "resource"
        source.queue_name = "shared"
        source.core_count = 24
        source.node_count = 1
        source.walltime_limit = 30
        source.total_physical_memory = 1024
        source.save()
        models.RemoteExecution.objects.create(
            run=source,
            airavata_experiment_id="completed-epolyscat",
            airavata_experiment_status="COMPLETED",
        )
        request = RequestFactory().post(
            f"/epolyscat_django_app/api/runs/{source.id}/workflow_continuation/",
            data={},
            content_type="application/json",
        )
        request.user = user
        viewset = views.RunViewSet()
        viewset.request = request
        viewset.get_object = mock.Mock(return_value=source)

        response = viewset.workflow_continuation(request, pk=source.id)

        self.assertIn("workflow_parent_run_id", response.data)
        parent = models.Run.objects.get(pk=response.data["workflow_parent_run_id"])
        child_runs = list(parent.workflow_steps.order_by("workflow_step_order"))
        self.assertEqual(parent.workflow_source_run, source)
        self.assertEqual(parent.run_mode, "workflow")
        self.assertTrue(parent.workflow_metadata["continuation"])
        self.assertEqual(
            parent.workflow_metadata["importedSource"],
            {
                "runId": source.id,
                "stage": "ePolyScat_Run",
                "application": "ePolyScat",
            },
        )
        self.assertEqual(
            [(child.workflow_step_order, child.workflow_stage) for child in child_runs],
            [(3, "Analysis")],
        )
        self.assertEqual(response.data["next_child_run_id"], child_runs[0].id)
        self.assertEqual(response.data["source_run_id"], source.id)
        self.assertIsNone(source.parent_run_id)

    def test_workflow_continuation_post_rejects_incomplete_run(self):
        user = get_user_model().objects.create_user(username="continue-incomplete")
        source = self.create_run(user)
        source.module_application = "ePolyScat"
        source.save(update_fields=["module_application"])
        models.RemoteExecution.objects.create(
            run=source,
            airavata_experiment_id="failed-epolyscat",
            airavata_experiment_status="FAILED",
        )
        request = RequestFactory().post(
            f"/epolyscat_django_app/api/runs/{source.id}/workflow_continuation/",
            data={},
            content_type="application/json",
        )
        request.user = user
        viewset = views.RunViewSet()
        viewset.request = request
        viewset.get_object = mock.Mock(return_value=source)

        with self.assertRaises(exceptions.ValidationError):
            viewset.workflow_continuation(request, pk=source.id)

        self.assertFalse(models.Run.objects.filter(workflow_source_run=source).exists())

    def test_workflow_continuation_presentation_includes_imported_source_stage(self):
        user = get_user_model().objects.create_user(username="continue-presentation")
        source = self.create_run(user)
        source.module_application = "ePolyScat"
        source.save(update_fields=["module_application"])
        parent = models.Run.objects.create(
            name="Continuation of source",
            owner=user,
            run_mode="workflow",
            workflow_source_run=source,
            workflow_metadata={
                "isWorkflowPlan": True,
                "continuation": True,
                "importedSource": {
                    "runId": source.id,
                    "stage": "ePolyScat_Run",
                    "application": "ePolyScat",
                },
            },
        )
        models.Run.objects.create(
            name="Analysis child",
            owner=user,
            run_mode="workflow",
            workflow_stage="Analysis",
            utility_application="CnvMath",
            parent_run=parent,
            workflow_step_order=3,
            workflow_step_status="pending",
        )
        request = RequestFactory().get(
            f"/epolyscat_django_app/api/runs/{parent.id}/"
        )
        request.user = user

        data = serializers.RunSerializer(parent, context={"request": request}).data
        stages = data["presentation"]["stages"]

        self.assertEqual(
            [(stage["id"], stage["state"], stage["status"]) for stage in stages],
            [
                ("Data_Gen", "not_included", "not_included"),
                ("ePolyScat_Run", "complete", "imported"),
                ("Analysis", "pending", "pending"),
                ("Visualization", "pending", "pending"),
            ],
        )
        self.assertEqual(stages[1]["child_run_id"], source.id)

    def test_workflow_presentation_splits_post_processing_and_visualization(self):
        serializer = serializers.RunSerializer()

        stages = serializer._workflow_stages("Analysis")

        self.assertEqual(
            [(stage["id"], stage["label"]) for stage in stages],
            [
                ("Data_Gen", "Data Generation"),
                ("ePolyScat_Run", "ePolyScat Run"),
                ("Analysis", "Post-processing"),
                ("Visualization", "Visualization"),
            ],
        )
        self.assertEqual(stages[2]["state"], "active")
        self.assertEqual(stages[3]["state"], "pending")
        self.assertTrue(stages[3]["local_only"])

    @override_settings(GATEWAY_DATA_STORE_REMOTE_API="https://amos-gateway.org/")
    @mock.patch("epolyscat_django_app.views.remoteapi.call")
    def test_create_run_user_dir_uses_remote_api_without_null_experiment_id(self, mock_call):
        mock_call.return_value.json.return_value = {"path": "EPOLYSCAT_Runs/Run_1"}
        request = RequestFactory().post("/epolyscat_django_app/api/runs/")
        run = models.Run(name="run", directory="EPOLYSCAT_Runs/Run_1")

        views._create_run_user_dir(request, run)

        _, kwargs = mock_call.call_args
        self.assertNotIn("data", kwargs)
        self.assertEqual(kwargs["path_params"], {"path": "EPOLYSCAT_Runs/Run_1"})

    @override_settings(GATEWAY_DATA_STORE_REMOTE_API="https://amos-gateway.org/")
    @mock.patch("epolyscat_django_app.views.remoteapi.call")
    def test_save_run_user_file_uses_remote_api_without_null_experiment_id(self, mock_call):
        mock_call.return_value.json.return_value = {
            "uploaded": {"productUri": "airavata-dp://saved"}
        }
        request = RequestFactory().post("/epolyscat_django_app/api/runs/")
        request.airavata_client = SimpleNamespace(
            getDataProduct=mock.Mock(
                return_value=SimpleNamespace(productUri="airavata-dp://saved")
            )
        )
        request.authz_token = object()
        run = models.Run(name="run", directory="EPOLYSCAT_Runs/Run_1")

        views._save_run_user_file(
            request,
            run,
            StringIO("contents"),
            name="input.dat",
            content_type="text/plain",
        )

        _, kwargs = mock_call.call_args
        self.assertNotIn("data", kwargs)
        self.assertEqual(kwargs["path_params"], {"path": "EPOLYSCAT_Runs/Run_1"})

    def test_save_file_rejects_storage_file_without_data_product_uri(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().post("/epolyscat_django_app/api/runs/")
        request.user = user
        run = self.create_run(user)
        input_obj = models.Input.objects.create(type="files", run=run, name="linp")

        viewset = views.RunViewSet()

        with self.assertRaises(exceptions.ValidationError):
            viewset._save_file(request, run, {"name": "linp"}, input_obj)

    def test_new_run_can_be_created_without_legacy_experiment_container(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )

        run = models.Run.objects.create(
            name="new run",
            owner=user,
            airavata_project_id="project-id",
        )

        self.assertIsNone(run.experiment)
        self.assertIsNone(run.root)
        self.assertEqual(run.number, "")
        self.assertEqual(run.filepath, "")

    def test_run_list_includes_user_workflow_parent_without_legacy_experiment(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/")
        request.user = user
        request.query_params = request.GET
        workflow_parent = models.Run.objects.create(
            name="workflow parent",
            owner=user,
            run_mode="workflow",
            workflow_application="OpenMolcas",
            workflow_metadata={"isWorkflowPlan": True},
        )
        tutorial_like_run = models.Run.objects.create(
            name="tutorial",
            owner=None,
            run_mode="workflow",
            workflow_metadata={"isWorkflowPlan": True},
        )
        viewset = views.RunViewSet()
        viewset.request = request
        viewset.action = "list"

        result = viewset.get_queryset()

        self.assertIn(workflow_parent, result)
        self.assertNotIn(tutorial_like_run, result)

    def test_run_and_view_lists_accept_gateway_user_objects(self):
        user = get_user_model().objects.create_user(
            username="gateway-user@example.com",
            password="password",
        )
        gateway_user = SimpleNamespace(
            username=user.username,
            email=user.email,
            first_name="",
            last_name="",
        )
        run = models.Run.objects.create(
            name="gateway run",
            owner=user,
            run_mode="module",
            module_application="ePolyScat",
        )
        view = models.View.objects.create(
            name="Gateway View",
            owner=user,
            type="default",
        )

        runs_request = RequestFactory().get("/epolyscat_django_app/api/runs/")
        runs_request.user = gateway_user
        runs_request.query_params = runs_request.GET
        runs_viewset = views.RunViewSet()
        runs_viewset.request = runs_request
        runs_viewset.action = "list"

        views_request = RequestFactory().get("/epolyscat_django_app/api/views/")
        views_request.user = gateway_user
        views_request.query_params = views_request.GET
        views_viewset = views.ViewsViewSet()
        views_viewset.request = views_request
        views_viewset.action = "list"

        self.assertIn(run, runs_viewset.get_queryset())
        self.assertIn(view, views_viewset.get_queryset())

    def test_show_viewable_route_allows_dotted_file_names(self):
        action = next(
            action
            for action in views.RunViewSet.get_extra_actions()
            if action.__name__ == "show_viewable"
        )

        self.assertEqual(action.url_path, r"viewables/(?P<filename>[^/]+)")

    @mock.patch("epolyscat_django_app.views.user_storage.open_file")
    @mock.patch("epolyscat_django_app.views.user_run_file_exists")
    def test_open_run_file_uses_data_product_uri_when_available(
        self, mock_exists, mock_open_file
    ):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/viewables/linp/")
        request.user = user
        run = self.create_run(user)
        mock_exists.return_value = "airavata-dp://input"
        mock_open_file.return_value = StringIO("contents")

        opened = views.open_run_file(request, run, "linp")

        self.assertEqual(opened.read(), "contents")
        mock_open_file.assert_called_once_with(
            request, data_product_uri="airavata-dp://input"
        )

    @mock.patch("epolyscat_django_app.views.user_storage.user_file_exists")
    def test_user_run_file_exists_checks_run_directory(self, mock_user_file_exists):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/viewables/linp/")
        request.user = user
        run = self.create_run(user)
        mock_user_file_exists.return_value = "airavata-dp://input"

        data_product_uri = views.user_run_file_exists(request, run, "linp")

        self.assertEqual(data_product_uri, "airavata-dp://input")
        mock_user_file_exists.assert_called_once_with(
            request, "EPOLYSCAT_Runs/Run_1/linp"
        )

    @override_settings(GATEWAY_ID="test-gateway")
    @mock.patch("epolyscat_django_app.views.experiment_util.launch")
    def test_create_remote_execution_includes_total_physical_memory(
        self, mock_launch
    ):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().post("/epolyscat_django_app/api/runs/1/submit/")
        request.user = user
        request.authz_token = object()
        request.airavata_client = SimpleNamespace(
            getApplicationInterface=mock.Mock(
                return_value=SimpleNamespace(
                    applicationInputs=[],
                    applicationOutputs=[],
                )
            ),
            createExperiment=mock.Mock(return_value="experiment-id"),
            getComputeResource=mock.Mock(return_value=SimpleNamespace(hostName="cluster")),
        )
        run = self.create_run(user)
        run.experiment.airavata_project_id = "project-id"
        run.experiment.save()
        run.group_resource_profile_id = "group-id"
        run.compute_resource_id = "compute-id"
        run.queue_name = "normal"
        run.core_count = 8
        run.node_count = 2
        run.walltime_limit = 120
        run.total_physical_memory = 4096

        views.RunViewSet()._create_remote_execution(
            request,
            run,
            {},
            "app-id",
            {},
            is_tutorial=False,
        )

        experiment_model = request.airavata_client.createExperiment.call_args.args[2]
        scheduling = experiment_model.userConfigurationData.computationalResourceScheduling
        self.assertEqual(scheduling.totalPhysicalMemory, 4096)
        mock_launch.assert_called_once_with(request, "experiment-id")

    @override_settings(GATEWAY_ID="test-gateway")
    @mock.patch("epolyscat_django_app.views.experiment_util.launch")
    def test_gaussian_module_adds_staged_input_to_command_line(self, mock_launch):
        user = get_user_model().objects.create_user(username="gaussian-module")
        request = RequestFactory().post("/epolyscat_django_app/api/runs/1/submit/")
        request.user = user
        request.authz_token = object()
        application_inputs = [
            SimpleNamespace(
                name="Calculation_Type",
                type=0,
                value=None,
                isRequired=True,
                requiredToAddedToCommandLine=True,
            ),
            SimpleNamespace(
                name="EPOLYSCAT_Application_Module",
                type=0,
                value=None,
                isRequired=True,
                requiredToAddedToCommandLine=True,
            ),
            SimpleNamespace(
                name="Gaussian_Inputs",
                type=4,
                value=None,
                isRequired=True,
                requiredToAddedToCommandLine=False,
            ),
        ]
        request.airavata_client = SimpleNamespace(
            getApplicationInterface=mock.Mock(
                return_value=SimpleNamespace(
                    applicationInputs=application_inputs,
                    applicationOutputs=[],
                )
            ),
            createExperiment=mock.Mock(return_value="experiment-id"),
            getComputeResource=mock.Mock(
                return_value=SimpleNamespace(hostName="cluster")
            ),
        )
        run = self.create_run(user)
        run.run_mode = "module"
        run.module_application = "Gaussian16"
        run.experiment.airavata_project_id = "project-id"
        run.experiment.save()
        run.group_resource_profile_id = "group-id"
        run.compute_resource_id = "compute-id"

        views.RunViewSet()._create_remote_execution(
            request,
            run,
            {},
            "app-id",
            {
                "Calculation_Type": "MODULE",
                "EPOLYSCAT_Application_Module": "Gaussian16",
                "Gaussian_Inputs": "airavata-dp://gaussian-input",
            },
            is_tutorial=False,
        )

        experiment = request.airavata_client.createExperiment.call_args.args[2]
        gaussian_input = next(
            inp for inp in experiment.experimentInputs
            if inp.name == "Gaussian_Inputs"
        )
        self.assertTrue(gaussian_input.requiredToAddedToCommandLine)
        mock_launch.assert_called_once_with(request, "experiment-id")

    @override_settings(GATEWAY_ID="test-gateway")
    @mock.patch("epolyscat_django_app.views.experiment_util.launch")
    def test_gaussian_workflow_adds_application_and_input_to_command_line(
        self, mock_launch
    ):
        user = get_user_model().objects.create_user(username="gaussian-workflow")
        request = RequestFactory().post("/epolyscat_django_app/api/runs/1/submit/")
        request.user = user
        request.authz_token = object()
        application_inputs = [
            SimpleNamespace(
                name="Calculation_Type",
                type=0,
                value=None,
                isRequired=True,
                requiredToAddedToCommandLine=True,
            ),
            SimpleNamespace(
                name="Application_Workflow",
                type=0,
                value=None,
                isRequired=True,
                requiredToAddedToCommandLine=True,
            ),
            SimpleNamespace(
                name="Data_Gen",
                type=0,
                value=None,
                isRequired=True,
                requiredToAddedToCommandLine=False,
            ),
            SimpleNamespace(
                name="Gaussian_Inputs",
                type=4,
                value=None,
                isRequired=True,
                requiredToAddedToCommandLine=False,
            ),
        ]
        request.airavata_client = SimpleNamespace(
            getApplicationInterface=mock.Mock(
                return_value=SimpleNamespace(
                    applicationInputs=application_inputs,
                    applicationOutputs=[],
                )
            ),
            createExperiment=mock.Mock(return_value="experiment-id"),
            getComputeResource=mock.Mock(
                return_value=SimpleNamespace(hostName="cluster")
            ),
        )
        run = self.create_run(user)
        run.run_mode = "workflow"
        run.workflow_stage = "Data_Gen"
        run.workflow_application = "Gaussian16"
        run.experiment.airavata_project_id = "project-id"
        run.experiment.save()
        run.group_resource_profile_id = "group-id"
        run.compute_resource_id = "compute-id"

        views.RunViewSet()._create_remote_execution(
            request,
            run,
            {},
            "app-id",
            {
                "Calculation_Type": "WORKFLOW",
                "Application_Workflow": "Data_Gen",
                "Data_Gen": "Gaussian16",
                "Gaussian_Inputs": "airavata-dp://gaussian-input",
            },
            is_tutorial=False,
        )

        experiment = request.airavata_client.createExperiment.call_args.args[2]
        command_line_flags = {
            inp.name: inp.requiredToAddedToCommandLine
            for inp in experiment.experimentInputs
        }
        self.assertTrue(command_line_flags["Data_Gen"])
        self.assertTrue(command_line_flags["Gaussian_Inputs"])
        mock_launch.assert_called_once_with(request, "experiment-id")

    @override_settings(GATEWAY_ID="test-gateway")
    @mock.patch("epolyscat_django_app.views.experiment_util.launch")
    def test_direct_epolyscat_deployment_only_passes_command_input_file(
        self, mock_launch
    ):
        user = get_user_model().objects.create_user(username="direct-epolyscat")
        request = RequestFactory().post("/epolyscat_django_app/api/runs/1/submit/")
        request.user = user
        request.authz_token = object()
        application_inputs = [
            SimpleNamespace(
                name="Calculation_Type",
                type=0,
                value=None,
                isRequired=True,
                requiredToAddedToCommandLine=True,
            ),
            SimpleNamespace(
                name="EPOLYSCAT_Application_Module",
                type=0,
                value=None,
                isRequired=True,
                requiredToAddedToCommandLine=True,
            ),
            SimpleNamespace(
                name="ePolyScat_Input_Data",
                type=4,
                value=None,
                isRequired=True,
                requiredToAddedToCommandLine=True,
            ),
            SimpleNamespace(
                name="ePolyscat_Input_File",
                type=3,
                value=None,
                isRequired=True,
                requiredToAddedToCommandLine=True,
            ),
        ]
        request.airavata_client = SimpleNamespace(
            getApplicationInterface=mock.Mock(
                return_value=SimpleNamespace(
                    applicationInputs=application_inputs,
                    applicationOutputs=[],
                )
            ),
            createExperiment=mock.Mock(return_value="experiment-id"),
            getComputeResource=mock.Mock(
                return_value=SimpleNamespace(hostName="frontera")
            ),
        )
        run = self.create_run(user)
        run.run_mode = "module"
        run.module_application = "ePolyScat"
        run.experiment.airavata_project_id = "project-id"
        run.experiment.save()
        run.group_resource_profile_id = "group-id"
        run.compute_resource_id = "frontera"

        input_values = {
            "Calculation_Type": "MODULE",
            "EPOLYSCAT_Application_Module": "ePolyScat",
            "ePolyScat_Input_Data": "airavata-dp://molden",
            "ePolyscat_Input_File": "airavata-dp://commands",
        }
        views.RunViewSet()._create_remote_execution(
            request,
            run,
            input_values,
            "app-id",
            input_values,
            is_tutorial=False,
            deployment_executable_path="/opt/epolyscat/bin/ePolyScat",
        )

        experiment = request.airavata_client.createExperiment.call_args.args[2]
        command_line_flags = {
            inp.name: inp.requiredToAddedToCommandLine
            for inp in experiment.experimentInputs
        }
        self.assertEqual(
            command_line_flags,
            {
                "Calculation_Type": False,
                "EPOLYSCAT_Application_Module": False,
                "ePolyScat_Input_Data": False,
                "ePolyscat_Input_File": True,
            },
        )
        mock_launch.assert_called_once_with(request, "experiment-id")

    @override_settings(
        GATEWAY_ID="test-gateway",
        EPOLYSCAT={
            "EPOLYSCAT_APPLICATION_ID": "epolyscat-module",
            "GAUSSIAN16_APPLICATION_ID": "gaussian-module",
        },
    )
    @mock.patch("epolyscat_django_app.views.user_storage.save_input_file")
    @mock.patch("epolyscat_django_app.views.user_storage.open_file")
    def test_gaussian_module_uses_dedicated_interface_and_input_name(
        self,
        mock_open_file,
        mock_save_input_file,
    ):
        user = get_user_model().objects.create_user(username="gaussian-route")
        request = RequestFactory().post("/epolyscat_django_app/api/runs/1/submit/")
        request.user = user
        request.authz_token = object()
        request.airavata_client = SimpleNamespace(
            getAllApplicationInterfaces=mock.Mock(
                return_value=[
                    SimpleNamespace(
                        applicationInterfaceId="epolyscat-interface",
                        applicationModules=["epolyscat-module"],
                    ),
                    SimpleNamespace(
                        applicationInterfaceId="gaussian-interface",
                        applicationModules=["gaussian-module"],
                    ),
                ]
            )
        )
        mock_open_file.return_value = StringIO("gaussian input")
        mock_save_input_file.return_value = SimpleNamespace(
            productUri="airavata-dp://staged-gaussian"
        )
        run = self.create_run(user)
        run.run_mode = "module"
        run.module_application = "Gaussian16"
        run.save()
        input_model = models.Input.objects.create(
            run=run,
            type="files",
            name="Gaussian_Input",
        )
        models.File.objects.create(
            input=input_model,
            name="gaussian.in",
            data_product_uri="airavata-dp://source-gaussian",
        )
        viewset = views.RunViewSet()
        viewset._create_remote_execution = mock.Mock()

        viewset._submit_single_run(request, run, is_tutorial=False)

        expected_inputs = {"Input_File": "airavata-dp://staged-gaussian"}
        viewset._create_remote_execution.assert_called_once_with(
            request,
            run,
            expected_inputs,
            "gaussian-interface",
            expected_inputs,
            False,
        )

    @override_settings(
        EPOLYSCAT={
            "EPOLYSCAT_APPLICATION_ID": "epolyscat-module",
            "GAUSSIAN16_APPLICATION_ID": "gaussian-module",
            "OPENMOLCAS_APPLICATION_ID": "openmolcas-module",
        }
    )
    def test_api_settings_exposes_gaussian_application_module(self):
        user = get_user_model().objects.create_user(username="api-settings")
        request = RequestFactory().get("/api/settings/")
        request.user = user

        response = views.api_settings(request)

        self.assertEqual(
            response.data["EPOLYSCAT"],
            {
                "EPOLYSCAT_APPLICATION_ID": "epolyscat-module",
                "GAUSSIAN16_APPLICATION_ID": "gaussian-module",
                "OPENMOLCAS_APPLICATION_ID": "openmolcas-module",
            },
        )

    @override_settings(
        GATEWAY_ID="test-gateway",
        EPOLYSCAT={
            "EPOLYSCAT_APPLICATION_ID": "epolyscat-module",
            "OPENMOLCAS_APPLICATION_ID": "openmolcas-module",
        },
    )
    @mock.patch("epolyscat_django_app.views.user_storage.save_input_file")
    @mock.patch("epolyscat_django_app.views.user_storage.open_file")
    def test_openmolcas_module_uses_dedicated_interface_and_input_name(
        self,
        mock_open_file,
        mock_save_input_file,
    ):
        user = get_user_model().objects.create_user(username="openmolcas-route")
        request = RequestFactory().post("/epolyscat_django_app/api/runs/1/submit/")
        request.user = user
        request.authz_token = object()
        request.airavata_client = SimpleNamespace(
            getAllApplicationInterfaces=mock.Mock(
                return_value=[
                    SimpleNamespace(
                        applicationInterfaceId="epolyscat-interface",
                        applicationModules=["epolyscat-module"],
                    ),
                    SimpleNamespace(
                        applicationInterfaceId="openmolcas-interface",
                        applicationModules=["openmolcas-module"],
                    ),
                ]
            )
        )
        mock_open_file.return_value = StringIO("openmolcas input")
        mock_save_input_file.return_value = SimpleNamespace(
            productUri="airavata-dp://staged-openmolcas"
        )
        run = self.create_run(user)
        run.run_mode = "module"
        run.module_application = "OpenMolcas"
        run.save()
        input_model = models.Input.objects.create(
            run=run,
            type="files",
            name="Molcas_Input",
        )
        models.File.objects.create(
            input=input_model,
            name="molcas.input",
            data_product_uri="airavata-dp://source-openmolcas",
        )
        viewset = views.RunViewSet()
        viewset._create_remote_execution = mock.Mock()

        viewset._submit_single_run(request, run, is_tutorial=False)

        expected_inputs = {
            "OpenMolcas input file": "airavata-dp://staged-openmolcas"
        }
        viewset._create_remote_execution.assert_called_once_with(
            request,
            run,
            expected_inputs,
            "openmolcas-interface",
            expected_inputs,
            False,
        )

    @override_settings(
        GATEWAY_ID="test-gateway",
        EPOLYSCAT={"EPOLYSCAT_APPLICATION_ID": "epolyscat-module"},
    )
    @mock.patch("epolyscat_django_app.views.user_storage.save_input_file")
    @mock.patch("epolyscat_django_app.views.user_storage.open_file")
    def test_epolyscat_module_passes_selected_deployment_entrypoint_to_execution(
        self,
        mock_open_file,
        mock_save_input_file,
    ):
        user = get_user_model().objects.create_user(username="epolyscat-route")
        request = RequestFactory().post("/epolyscat_django_app/api/runs/1/submit/")
        request.user = user
        request.authz_token = object()
        request.airavata_client = SimpleNamespace(
            getAllApplicationInterfaces=mock.Mock(
                return_value=[
                    SimpleNamespace(
                        applicationInterfaceId="epolyscat-interface",
                        applicationModules=["epolyscat-module"],
                    ),
                ]
            ),
            getAllApplicationDeployments=mock.Mock(
                return_value=[
                    SimpleNamespace(
                        appModuleId="epolyscat-module",
                        computeHostId="frontera",
                        executablePath="/opt/epolyscat/bin/ePolyScat",
                    ),
                ]
            ),
        )
        mock_open_file.return_value = StringIO("input")
        mock_save_input_file.side_effect = [
            SimpleNamespace(productUri="airavata-dp://staged-molden"),
            SimpleNamespace(productUri="airavata-dp://staged-commands"),
        ]
        run = self.create_run(user)
        run.run_mode = "module"
        run.module_application = "ePolyScat"
        run.compute_resource_id = "frontera"
        run.save()
        data_input = models.Input.objects.create(
            run=run,
            type="files",
            name="ePolyScat_Input_Data",
        )
        models.File.objects.create(
            input=data_input,
            name="water.scf.molden",
            data_product_uri="airavata-dp://source-molden",
        )
        command_input = models.Input.objects.create(
            run=run,
            type="files",
            name="ePolyscat_Input_File",
        )
        models.File.objects.create(
            input=command_input,
            name="test03.inp",
            data_product_uri="airavata-dp://source-commands",
        )
        viewset = views.RunViewSet()
        viewset._create_remote_execution = mock.Mock()

        viewset._submit_single_run(request, run, is_tutorial=False)

        self.assertEqual(
            viewset._create_remote_execution.call_args.kwargs,
            {"deployment_executable_path": "/opt/epolyscat/bin/ePolyScat"},
        )
        positional = viewset._create_remote_execution.call_args.args
        self.assertEqual(positional[3], "epolyscat-interface")
        self.assertEqual(
            positional[4]["ePolyscat_Input_File"],
            "airavata-dp://staged-commands",
        )

    @mock.patch("epolyscat_django_app.views.user_storage.save_input_file")
    @mock.patch("epolyscat_django_app.views.user_storage.open_file")
    def test_resubmit_calls_create_remote_execution_with_current_signature(
        self,
        mock_open_file,
        mock_save_input_file,
    ):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().post("/epolyscat_django_app/api/runs/1/resubmit/")
        request.user = user
        request.data = {}
        run = self.create_run(user)
        run.inpc_data_product_uri = "airavata-dp://inpc"
        run.save()
        models.RemoteExecution.objects.create(
            run=run,
            airavata_experiment_id="old-experiment-id",
            job_id="job-1",
        )
        mock_open_file.return_value = StringIO("inpc")
        mock_save_input_file.return_value = SimpleNamespace(
            productUri="airavata-dp://saved-inpc"
        )
        input_serializer = mock.Mock()
        input_serializer.is_valid.return_value = None
        input_serializer.save.return_value = run
        response_serializer = mock.Mock()
        response_serializer.data = {"id": run.id}
        viewset = views.RunViewSet()
        viewset.get_object = mock.Mock(return_value=run)
        viewset.get_serializer = mock.Mock(
            side_effect=[input_serializer, response_serializer]
        )
        viewset._get_eployscat_app_interface_id = mock.Mock(return_value="app-id")
        viewset._create_remote_execution = mock.Mock()

        result = viewset.resubmit(request, pk=run.id)

        input_values = {
            "Input-File": "airavata-dp://saved-inpc",
            "Previous_JobID_Restart": "job-1",
        }
        viewset._create_remote_execution.assert_called_once_with(
            request,
            run,
            input_values,
            "app-id",
            input_values,
            False,
        )
        self.assertEqual(result.data, {"id": run.id})

    @mock.patch("epolyscat_django_app.views.user_storage.list_experiment_dir")
    def test_get_output_files_merges_root_and_archive_outputs(self, mock_list_experiment_dir):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/get_output_files/")
        request.user = user
        run = self.create_run(user)
        models.RemoteExecution.objects.create(
            run=run,
            airavata_experiment_id="experiment-id",
        )
        root_files = [
            {"name": "job_1683512072.slurm", "data-product-uri": "airavata-dp://slurm"},
            {"name": "ePolyScat.stdout", "data-product-uri": "airavata-dp://stdout"},
        ]
        archive_files = [
            {"name": "gaussian.log", "data-product-uri": "airavata-dp://gaussian-log"},
            {"name": "molden.dat", "data-product-uri": "airavata-dp://molden"},
        ]

        def list_experiment_dir(_request, _experiment_id, path=None):
            return ([], archive_files if path == "ARCHIVE" else root_files)

        mock_list_experiment_dir.side_effect = list_experiment_dir
        viewset = views.RunViewSet()
        viewset.get_object = mock.Mock(return_value=run)

        response = viewset.get_output_files(request, pk=run.id)

        self.assertEqual(
            [file_data["name"] for file_data in response.data],
            ["job_1683512072.slurm", "ePolyScat.stdout", "gaussian.log", "molden.dat"],
        )
        self.assertEqual(mock_list_experiment_dir.call_count, 2)

    @mock.patch("epolyscat_django_app.views.user_storage.list_experiment_dir")
    def test_get_output_files_recurses_experiment_output_directories(self, mock_list_experiment_dir):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/get_output_files/")
        request.user = user
        run = self.create_run(user)
        models.RemoteExecution.objects.create(
            run=run,
            airavata_experiment_id="experiment-id",
        )

        def list_experiment_dir(_request, _experiment_id, path=""):
            tree = {
                "": (
                    [{"name": "ARCHIVE"}, {"name": "scratch"}],
                    [{"name": "job_1683512072.slurm", "data-product-uri": "airavata-dp://slurm"}],
                ),
                "ARCHIVE": (
                    [{"name": "Gaussian"}],
                    [{"name": "ePolyScat.stdout", "data-product-uri": "airavata-dp://stdout"}],
                ),
                "ARCHIVE/Gaussian": (
                    [],
                    [
                        {"name": "gaussian.log", "data-product-uri": "airavata-dp://gaussian-log"},
                        {"name": "molden.dat", "data-product-uri": "airavata-dp://molden"},
                    ],
                ),
                "scratch": (
                    [],
                    [{"name": "gaussian.inp", "data-product-uri": "airavata-dp://gaussian-inp"}],
                ),
            }
            return tree[path]

        mock_list_experiment_dir.side_effect = list_experiment_dir
        viewset = views.RunViewSet()
        viewset.get_object = mock.Mock(return_value=run)

        response = viewset.get_output_files(request, pk=run.id)

        self.assertEqual(
            [file_data["name"] for file_data in response.data],
            [
                "job_1683512072.slurm",
                "ePolyScat.stdout",
                "gaussian.log",
                "molden.dat",
                "gaussian.inp",
            ],
        )

    @mock.patch("epolyscat_django_app.views.user_storage.list_experiment_dir")
    def test_workflow_output_binding_returns_backend_selected_file(self, mock_list_experiment_dir):
        user = get_user_model().objects.create_user(username="binding-user")
        run = self.create_run(user)
        run.module_application = "Gaussian16"
        run.save(update_fields=["module_application"])
        models.RemoteExecution.objects.create(
            run=run,
            airavata_experiment_id="gaussian-experiment",
        )
        mock_list_experiment_dir.return_value = (
            [],
            [
                {"name": "job.slurm", "data-product-uri": "airavata-dp://job"},
                {"name": "gaussian.log", "data-product-uri": "airavata-dp://log"},
            ],
        )
        request = RequestFactory().get(
            f"/epolyscat_django_app/api/runs/{run.id}/workflow_output_binding/",
            {"targetStageId": "ePolyScat_Run"},
        )
        request.user = user
        viewset = views.RunViewSet()
        viewset.request = request
        viewset.get_object = mock.Mock(return_value=run)

        response = viewset.workflow_output_binding(request, pk=run.id)

        self.assertEqual(response.data["status"], "ready")
        self.assertEqual(response.data["source_application"], "Gaussian16")
        self.assertEqual(response.data["selected"]["name"], "gaussian.log")
        self.assertEqual(response.data["input_file_name"], "ePolyScat_Input_Data")

    @override_settings(
        GATEWAY_ID="gateway",
        EPOLYSCAT={"EPOLYSCAT_APPLICATION_ID": "epolyscat-module"},
    )
    def test_runtime_audit_action_uses_airavata_runtime_configuration(self):
        request = RequestFactory().get(
            "/epolyscat_django_app/api/runs/runtime_audit/"
        )
        request.user = get_user_model().objects.create_user(username="audit-user")
        request.authz_token = "token"
        request.airavata_client = SimpleNamespace(
            getAllAppModules=mock.Mock(return_value=[]),
            getAllApplicationInterfaces=mock.Mock(return_value=[]),
            getAllApplicationDeployments=mock.Mock(return_value=[]),
        )
        viewset = views.RunViewSet()
        viewset.request = request

        response = viewset.runtime_audit(request)

        self.assertFalse(response.data["ready"])
        self.assertEqual(
            [utility["id"] for utility in response.data["utilities"]],
            [
                "CnvMath",
                "CnvMatLab",
                "CnvLinFull",
                "MoldenMerge",
                "NRFPAD",
                "Cube2igor",
            ],
        )
        request.airavata_client.getAllAppModules.assert_called_once_with(
            "token", "gateway"
        )

    def test_submit_normalizes_workflow_data_generation_inputs_for_airavata(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().post("/epolyscat_django_app/api/runs/1/submit/")
        request.user = user
        request.data = {}
        run = self.create_run(user)
        run.run_mode = "workflow"
        run.workflow_stage = "data-generation"
        run.workflow_application = "OpenMolcas"
        run.save()
        models.Input.objects.create(
            run=run,
            type="radio input",
            name="Calculation_Type",
            value="WORKFLOW",
        )
        models.Input.objects.create(
            run=run,
            type="radio input",
            name="Application_Workflow",
            value="data-generation",
        )
        models.Input.objects.create(
            run=run,
            type="parameter",
            name="Workflow_Application",
            value="OpenMolcas",
        )
        molcas_input = models.Input.objects.create(
            run=run,
            type="files",
            name="Molcas_Input",
        )
        models.File.objects.create(
            input=molcas_input,
            name="molcas.input",
            data_product_uri="airavata-dp://source-openmolcas",
        )
        serializer = mock.Mock()
        serializer.is_valid.return_value = None
        serializer.save.return_value = run
        serializer.data = {"is_tutorial": False}
        response_serializer = mock.Mock()
        response_serializer.data = {"id": run.id}
        viewset = views.RunViewSet()
        viewset.get_object = mock.Mock(return_value=run)
        viewset.get_serializer = mock.Mock(side_effect=[serializer, response_serializer])
        viewset._get_run_app_interface_id = mock.Mock(return_value="app-id")
        viewset._create_remote_execution = mock.Mock()

        with mock.patch(
            "epolyscat_django_app.views.user_storage.open_file",
            return_value=StringIO("openmolcas input"),
        ), mock.patch(
            "epolyscat_django_app.views.user_storage.save_input_file",
            return_value=SimpleNamespace(
                productUri="airavata-dp://staged-openmolcas"
            ),
        ):
            result = viewset.submit(request, pk=run.id)

        expected_inputs = {
            "OpenMolcas input file": "airavata-dp://staged-openmolcas",
        }
        viewset._create_remote_execution.assert_called_once_with(
            request,
            run,
            expected_inputs,
            "app-id",
            expected_inputs,
            False,
        )
        self.assertEqual(result.data, {"id": run.id})

    def test_airavata_input_builder_supports_module_workflow_and_utility(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        run = self.create_run(user)

        run.run_mode = "module"
        run.module_application = "ePolyScat"
        module_inputs = views.build_airavata_input_values(
            run,
            {
                "ePolyScat_Input_Data": "airavata-dp://data",
                "ePolyscat_Input_File": "airavata-dp://commands",
                "molden.dat": "airavata-dp://molden",
            },
        )
        self.assertEqual(module_inputs["Calculation_Type"], "MODULE")
        self.assertEqual(module_inputs["EPOLYSCAT_Application_Module"], "ePolyScat")
        self.assertEqual(module_inputs["ePolyScat_Input_Data"], "airavata-dp://data")
        self.assertEqual(module_inputs["ePolyscat_Input_File"], "airavata-dp://commands")
        self.assertEqual(module_inputs["molden.dat"], "airavata-dp://molden")

        run.run_mode = "workflow"
        run.workflow_stage = "ePolyScat_Run"
        workflow_inputs = views.build_airavata_input_values(
            run,
            {
                "ePolyScat_Input_Data": "airavata-dp://data",
                "ePolyscat_Input_File": "airavata-dp://commands",
            },
        )
        self.assertEqual(workflow_inputs["Calculation_Type"], "WORKFLOW")
        self.assertEqual(workflow_inputs["Application_Workflow"], "ePolyScat_Run")
        self.assertEqual(workflow_inputs["ePolyScat_Input_Data"], "airavata-dp://data")
        self.assertNotIn("EPOLYSCAT_Application_Module", workflow_inputs)

        run.run_mode = "utility"
        run.utility_application = "MoldenMerge"
        utility_inputs = views.build_airavata_input_values(
            run,
            {"molden.dat": "airavata-dp://molden"},
        )
        self.assertEqual(utility_inputs["Calculation_Type"], "UTILITY")
        self.assertEqual(utility_inputs["Application_Utility"], "MoldenMerge")
        self.assertEqual(utility_inputs["molden.dat"], "airavata-dp://molden")

        run.utility_application = "CnvLinFull"
        utility_inputs = views.build_airavata_input_values(
            run,
            {
                "DumpOut": "airavata-dp://dump",
                "ePolyscat_Input_Data": "airavata-dp://existing",
            },
        )
        self.assertEqual(utility_inputs["Application_Utility"], "CnvLinFull")
        self.assertEqual(
            utility_inputs["ePolyscat_Input_Data"],
            "airavata-dp://existing,airavata-dp://dump",
        )

        run.run_mode = "workflow"
        run.workflow_stage = "Analysis"
        run.utility_application = "NRFPAD"
        workflow_analysis_inputs = views.build_airavata_input_values(
            run,
            {"Cross_Section_Input_File": "airavata-dp://orient"},
        )
        self.assertEqual(workflow_analysis_inputs["Calculation_Type"], "WORKFLOW")
        self.assertEqual(workflow_analysis_inputs["Application_Workflow"], "Analysis")
        self.assertEqual(workflow_analysis_inputs["Application_Utility"], "NRFPAD")
        self.assertEqual(
            workflow_analysis_inputs["ePolyscat_Input_Data"],
            "airavata-dp://orient",
        )

    def test_run_serializer_exposes_workflow_presentation_metadata(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/")
        request.user = user
        run = self.create_run(user)
        run.run_mode = "workflow"
        run.workflow_stage = "bound"
        run.workflow_application = "OpenMolcas"
        run.workflow_metadata = {"selected_file": "target"}
        run.save()

        data = serializers.RunSerializer(run, context={"request": request}).data

        self.assertEqual(data["run_mode"], "workflow")
        self.assertEqual(data["workflow_stage"], "bound")
        self.assertEqual(data["workflow_application"], "OpenMolcas")
        self.assertEqual(data["workflow_metadata"], {"selected_file": "target"})
        self.assertEqual(data["presentation"]["mode"], "workflow")
        self.assertEqual(data["presentation"]["subtitle"], "Workflow/Bound")
        self.assertEqual(data["presentation"]["active_stage_id"], "bound")
        self.assertEqual(
            [application["id"] for application in data["presentation"]["applications"]],
            ["Gaussian16", "OpenMolcas"],
        )

    def test_run_presentation_action_returns_module_read_model(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().get(
            "/epolyscat_django_app/api/runs/1/presentation/"
        )
        request.user = user
        run = self.create_run(user)

        viewset = views.RunViewSet()
        viewset.get_object = mock.Mock(return_value=run)
        viewset.get_serializer = mock.Mock(
            return_value=serializers.RunSerializer(run, context={"request": request})
        )

        result = viewset.presentation(request, pk=run.id)

        self.assertEqual(result.data["mode"], "module")
        self.assertEqual(result.data["subtitle"], "Modules/EPOLYSCAT_DMAT")
        self.assertIn("target_states", result.data)
        self.assertEqual(
            result.data["plottable_file_names"],
            ["spec_total", "spec_partial", "expec", "Laser", "eig"],
        )

    def test_run_presentation_files_follow_module_application(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/")
        request.user = user
        run = self.create_run(user)
        run.run_mode = "module"
        run.module_application = "Gaussian16"
        run.save()

        data = serializers.RunSerializer(run, context={"request": request}).data

        presentation = data["presentation"]
        flattened_inputs = [
            filename
            for column in presentation["target_states"]["input_columns"]
            for filename in column
        ]
        self.assertEqual(presentation["subtitle"], "Modules/Gaussian16")
        self.assertEqual(presentation["selected_file"], "Gaussian_Input")
        self.assertEqual(flattened_inputs, ["Gaussian_Input"])
        self.assertIn("gaussian.log", presentation["target_states"]["output_files"])
        self.assertNotIn("ns_001.c", flattened_inputs)
        self.assertNotIn("ePolyScat_dmat.log", presentation["target_states"]["output_files"])

    def test_workflow_parent_presentation_files_follow_active_child_stage(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/")
        request.user = user
        parent = self.create_run(user)
        parent.run_mode = "workflow"
        parent.workflow_application = "OpenMolcas"
        parent.save()
        child_runs = views.RunViewSet()._ensure_workflow_child_runs(parent)
        child_runs[0].workflow_step_status = "submitted"
        child_runs[0].save()

        data = serializers.RunSerializer(parent, context={"request": request}).data

        presentation = data["presentation"]
        flattened_inputs = [
            filename
            for column in presentation["target_states"]["input_columns"]
            for filename in column
        ]
        self.assertEqual(presentation["active_stage_id"], "Data_Gen")
        self.assertEqual(presentation["selected_file"], "Molcas_Input")
        self.assertEqual(flattened_inputs, ["Molcas_Input"])
        self.assertIn("molcas.log", presentation["target_states"]["output_files"])
        self.assertNotIn("ns_001.c", flattened_inputs)

    def test_workflow_parent_creates_ordered_child_runs(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        parent = self.create_run(user)
        parent.run_mode = "workflow"
        parent.workflow_application = "OpenMolcas"
        parent.workflow_metadata = {
            "analysisApplications": ["CnvMath", "CnvLinFull"],
        }
        parent.group_resource_profile_id = "group-id"
        parent.compute_resource_id = "compute-id"
        parent.queue_name = "normal"
        parent.core_count = 4
        parent.node_count = 1
        parent.walltime_limit = 60
        parent.total_physical_memory = 2048
        parent.save()
        parent_input = models.Input.objects.create(
            run=parent,
            type="files",
            name="Molcas_Input",
            value=None,
        )
        models.File.objects.create(
            input=parent_input,
            name="molcas.input",
            data_product_uri="airavata-dp://molcas",
        )

        child_runs = views.RunViewSet()._ensure_workflow_child_runs(parent)

        self.assertEqual(
            [
                (child.workflow_step_order, child.workflow_stage)
                for child in child_runs
            ],
            [
                (1, "Data_Gen"),
                (2, "ePolyScat_Run"),
                (3, "Analysis"),
                (4, "Analysis"),
            ],
        )
        self.assertEqual(child_runs[0].workflow_application, "OpenMolcas")
        self.assertEqual(child_runs[2].utility_application, "CnvMath")
        self.assertEqual(child_runs[3].utility_application, "CnvLinFull")
        self.assertEqual(child_runs[0].parent_run, parent)
        self.assertEqual(child_runs[0].group_resource_profile_id, "group-id")
        self.assertEqual(child_runs[0].inputs.get().files.get().data_product_uri, "airavata-dp://molcas")

    def test_multiple_analysis_children_render_as_one_post_processing_step(self):
        user = get_user_model().objects.create_user(
            username="multi-analysis@example.com",
            password="password",
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/")
        request.user = user
        parent = self.create_run(user)
        parent.run_mode = "workflow"
        parent.workflow_application = "OpenMolcas"
        parent.workflow_metadata = {
            "analysisApplications": ["CnvMath", "CnvLinFull"],
        }
        parent.save()
        child_runs = views.RunViewSet()._ensure_workflow_child_runs(parent)

        data = serializers.RunSerializer(parent, context={"request": request}).data

        stages = data["presentation"]["stages"]
        self.assertEqual(
            [stage["id"] for stage in stages],
            ["Data_Gen", "ePolyScat_Run", "Analysis", "Visualization"],
        )
        analysis = stages[2]
        analysis_children = [
            child for child in child_runs if child.workflow_stage == "Analysis"
        ]
        self.assertEqual(
            analysis["child_run_ids"],
            [child.id for child in analysis_children],
        )
        self.assertEqual(
            analysis["applications"],
            ["CnvMath", "CnvLinFull"],
        )

    def test_workflow_submit_delegates_to_first_child_run(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().post("/epolyscat_django_app/api/runs/1/submit/")
        request.user = user
        parent = self.create_run(user)
        parent.run_mode = "workflow"
        parent.workflow_application = "OpenMolcas"
        parent.save()
        viewset = views.RunViewSet()
        viewset._submit_single_run = mock.Mock()

        result = viewset._submit_workflow_run(request, parent, is_tutorial=False)

        first_child = parent.workflow_steps.get(workflow_step_order=1)
        viewset._submit_single_run.assert_called_once_with(
            request,
            first_child,
            is_tutorial=False,
        )
        parent.refresh_from_db()
        self.assertEqual(parent.workflow_metadata["active_child_run_id"], first_child.id)
        self.assertEqual(result, parent)

    def test_workflow_child_submit_requires_previous_steps_to_finish(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().post("/epolyscat_django_app/api/runs/2/submit/")
        request.user = user
        request.data = {}
        parent = self.create_run(user)
        parent.run_mode = "workflow"
        parent.workflow_application = "OpenMolcas"
        parent.save()
        child_runs = views.RunViewSet()._ensure_workflow_child_runs(parent)
        second_child = child_runs[1]
        viewset = views.RunViewSet()
        viewset.get_object = mock.Mock(return_value=second_child)
        serializer = mock.Mock()
        serializer.is_valid.return_value = None
        serializer.save.return_value = second_child
        serializer.data = {"is_tutorial": False}
        viewset.get_serializer = mock.Mock(return_value=serializer)
        viewset._submit_single_run = mock.Mock()

        with self.assertRaises(exceptions.ValidationError) as context:
            viewset.submit(request, pk=second_child.id)

        self.assertIn("Data Generation", str(context.exception))
        viewset._submit_single_run.assert_not_called()

    def test_workflow_presentation_includes_child_run_statuses(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/")
        request.user = user
        parent = self.create_run(user)
        parent.run_mode = "workflow"
        parent.workflow_application = "OpenMolcas"
        parent.save()
        child_runs = views.RunViewSet()._ensure_workflow_child_runs(parent)
        child_runs[0].workflow_step_status = "submitted"
        child_runs[0].save()

        data = serializers.RunSerializer(parent, context={"request": request}).data

        stages = data["presentation"]["stages"]
        self.assertEqual(stages[0]["child_run_id"], child_runs[0].id)
        self.assertEqual(stages[0]["application"], "OpenMolcas")
        self.assertEqual(stages[0]["state"], "active")
        self.assertEqual(stages[1]["state"], "pending")

    def test_workflow_presentation_marks_finished_child_step_complete(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/")
        request.user = user
        parent = self.create_run(user)
        parent.run_mode = "workflow"
        parent.workflow_application = "OpenMolcas"
        parent.save()
        child_runs = views.RunViewSet()._ensure_workflow_child_runs(parent)
        child_runs[0].workflow_step_status = "submitted"
        child_runs[0].save()
        models.RemoteExecution.objects.create(
            run=child_runs[0],
            airavata_experiment_id="exp-1",
            airavata_experiment_status="COMPLETED",
        )

        data = serializers.RunSerializer(parent, context={"request": request}).data

        stages = data["presentation"]["stages"]
        self.assertEqual(stages[0]["state"], "complete")
        self.assertEqual(stages[0]["status"], "complete")
        self.assertEqual(stages[1]["state"], "pending")

    def test_workflow_visualization_becomes_ready_after_post_processing_completes(self):
        user = get_user_model().objects.create_user(
            username="visualization-ready@example.com",
            password="password",
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/")
        request.user = user
        parent = self.create_run(user)
        parent.run_mode = "workflow"
        parent.workflow_application = "OpenMolcas"
        parent.save()
        child_runs = views.RunViewSet()._ensure_workflow_child_runs(parent)
        analysis_child = next(
            child for child in child_runs if child.workflow_stage == "Analysis"
        )
        analysis_child.workflow_step_status = "complete"
        analysis_child.save(update_fields=["workflow_step_status"])

        data = serializers.RunSerializer(parent, context={"request": request}).data

        visualization = data["presentation"]["stages"][-1]
        self.assertEqual(visualization["id"], "Visualization")
        self.assertEqual(visualization["state"], "active")
        self.assertEqual(visualization["status"], "ready")
        self.assertTrue(visualization["local_only"])
        self.assertEqual(
            visualization["source_child_run_id"],
            analysis_child.id,
        )

    def test_workflow_presentation_exposes_parent_progress_status(self):
        user = get_user_model().objects.create_user(
            username="user@example.com",
            password="password",
        )
        request = RequestFactory().get("/epolyscat_django_app/api/runs/1/")
        request.user = user
        parent = self.create_run(user)
        parent.run_mode = "workflow"
        parent.workflow_application = "OpenMolcas"
        parent.save()
        child_runs = views.RunViewSet()._ensure_workflow_child_runs(parent)
        child_runs[0].workflow_step_status = "submitted"
        child_runs[0].save()
        parent.workflow_metadata = {
            "workflow_state": "running",
            "active_child_run_id": child_runs[0].id,
        }
        parent.save()

        data = serializers.RunSerializer(parent, context={"request": request}).data

        self.assertEqual(data["presentation"]["workflow_state"], "running")
        self.assertEqual(data["presentation"]["active_child_run_id"], child_runs[0].id)
        self.assertEqual(data["presentation"]["active_child_label"], "Data Generation")
        self.assertEqual(data["presentation"]["active_child_status"], "submitted")
