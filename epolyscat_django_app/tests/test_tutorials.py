import json
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from django.contrib.auth import get_user_model
from django.db import IntegrityError, connection, transaction
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext
from rest_framework.test import APIRequestFactory, force_authenticate

from epolyscat_django_app import models, serializers, tutorials, views


class TutorialViewBackendTests(TestCase):
    def setUp(self):
        self.user = get_user_model().objects.create_user(
            username="tutorial-reader",
        )
        self.factory = APIRequestFactory()

    def authenticated_request(self, method, path, data=None):
        request = getattr(self.factory, method)(path, data=data, format="json")
        request.user = self.user
        force_authenticate(request, user=self.user)
        return request

    def test_loads_curated_manual_cases_in_manifest_order(self):
        self.assertTrue(
            hasattr(tutorials, "load_catalog"),
            "tutorials.load_catalog must load the packaged manual catalog",
        )

        cases = tutorials.load_catalog()

        self.assertEqual(
            [case["id"] for case in cases],
            ["test01", "test04", "test08", "test12", "test36"],
        )
        self.assertEqual(
            {case["convert_format"] for case in cases},
            {"molden", "gaussian", "gamess"},
        )

    def test_catalog_rejects_an_unapproved_case_set(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            catalog_root = Path(temporary_directory)
            manifest_path = catalog_root / "manifest.json"
            manifest_path.write_text(
                json.dumps({"version": 1, "cases": []}),
                encoding="utf-8",
            )

            with mock.patch.object(tutorials, "MANIFEST_PATH", manifest_path):
                with self.assertRaisesMessage(
                    ValueError,
                    "exactly the approved cases",
                ):
                    tutorials.load_catalog()

    def test_catalog_rejects_unsafe_artifact_paths(self):
        manifest = json.loads(tutorials.MANIFEST_PATH.read_text(encoding="utf-8"))
        manifest["cases"][0]["script"] = "../escape.inp"

        with tempfile.TemporaryDirectory() as temporary_directory:
            manifest_path = Path(temporary_directory) / "manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            with mock.patch.object(tutorials, "MANIFEST_PATH", manifest_path):
                with self.assertRaisesMessage(ValueError, "direct basenames"):
                    tutorials.load_catalog()

    def test_catalog_requires_existing_artifacts_and_supported_formats(self):
        original_manifest = json.loads(
            tutorials.MANIFEST_PATH.read_text(encoding="utf-8")
        )
        mutations = (
            ("runtime_data", "missing.dat", "declared artifact"),
            ("convert_format", "unsupported", "convert format"),
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            for field, value, message in mutations:
                with self.subTest(field=field):
                    manifest = json.loads(json.dumps(original_manifest))
                    manifest["cases"][0][field] = value
                    manifest_path = Path(temporary_directory) / f"{field}.json"
                    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

                    with mock.patch.object(tutorials, "MANIFEST_PATH", manifest_path):
                        with self.assertRaisesMessage(ValueError, message):
                            tutorials.load_catalog()

    def test_catalog_rejects_artifact_symlinks_outside_the_case_directory(self):
        manifest = json.loads(
            tutorials.MANIFEST_PATH.read_text(encoding="utf-8")
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            catalog_root = Path(temporary_directory) / "manual_tutorials"
            catalog_root.mkdir()
            manifest_path = catalog_root / "manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            external_artifact = Path(temporary_directory) / "external.inp"
            external_artifact.write_text("external", encoding="utf-8")

            for case in manifest["cases"]:
                case_path = catalog_root / case["id"]
                case_path.mkdir()
                for field in tutorials.ARTIFACT_FIELDS:
                    artifact_path = case_path / case[field]
                    if case["id"] == "test01" and field == "script":
                        artifact_path.symlink_to(external_artifact)
                    else:
                        artifact_path.write_text("fixture", encoding="utf-8")

            with mock.patch.object(
                tutorials,
                "MANUAL_TUTORIALS_DIR",
                catalog_root,
            ), mock.patch.object(tutorials, "MANIFEST_PATH", manifest_path):
                with self.assertRaisesMessage(ValueError, "inside its case"):
                    tutorials.load_catalog()

    def test_synchronizes_one_shared_view_and_five_cloneable_runs_idempotently(self):
        first_view = tutorials.ensure_tutorial_view()
        stale_run = models.Run.objects.create(
            owner=None,
            name="Old tutorial",
            number="old-tutorial",
            filepath="data/epolyscat/tutorials/old-tutorial",
        )
        first_view.runs.add(stale_run)
        second_view = tutorials.ensure_tutorial_view()
        extra_file_input = models.Input.objects.create(
            run=second_view.runs.order_by("pk").first(),
            type="files",
            name="Unexpected_Tutorial_Input",
        )
        models.File.objects.create(
            input=extra_file_input,
            name="unexpected.dat",
        )
        extra_parameter_input = models.Input.objects.create(
            run=second_view.runs.order_by("pk").first(),
            type="parameter",
            name="Unexpected_Tutorial_Parameter",
            value="must not reach clone hydration",
        )
        synchronized_view = tutorials.ensure_tutorial_view()

        self.assertEqual(synchronized_view.pk, first_view.pk)
        self.assertIsNone(first_view.owner)
        self.assertEqual(first_view.type, "tutorial")
        self.assertEqual(models.View.objects.filter(type="tutorial").count(), 1)

        runs = list(first_view.runs.order_by("name"))
        self.assertEqual(
            [run.name for run in runs],
            [
                "Test 01",
                "Test 04",
                "Test 08",
                "Test 12",
                "Test 36",
            ],
        )
        self.assertTrue(models.Run.objects.filter(pk=stale_run.pk).exists())
        self.assertFalse(synchronized_view.runs.filter(pk=stale_run.pk).exists())
        self.assertTrue(all(run.owner is None for run in runs))
        self.assertTrue(all(run.run_mode == "module" for run in runs))
        self.assertTrue(all(run.module_application == "ePolyScat" for run in runs))
        self.assertFalse(
            models.Input.objects.filter(pk=extra_file_input.pk).exists()
        )
        self.assertFalse(
            models.Input.objects.filter(pk=extra_parameter_input.pk).exists()
        )

        for run in runs:
            selector_values = {
                input_instance.name: input_instance.value
                for input_instance in run.inputs.exclude(type="files")
            }
            self.assertEqual(selector_values["Calculation_Type"], "MODULE")
            self.assertEqual(
                selector_values["EPOLYSCAT_Application_Module"],
                "ePolyScat",
            )
            case = next(
                catalog_case
                for catalog_case in tutorials.load_catalog()
                if catalog_case["id"] == run.workflow_metadata["tutorialCaseId"]
            )
            self.assertEqual(
                run.description,
                case["description"],
            )
            self.assertEqual(
                run.workflow_metadata["tutorialInputFormat"],
                case["input_format"],
            )
            self.assertEqual(
                run.workflow_metadata["tutorialConvertFormat"],
                case["convert_format"],
            )
            self.assertTrue(run.workflow_metadata["tutorialStandalone"])
            self.assertEqual(
                set(run.inputs.filter(type="files").values_list("name", flat=True)),
                {"ePolyscat_Input_File", "ePolyScat_Input_Data"},
            )
            self.assertEqual(
                list(
                    run.inputs.get(name="ePolyscat_Input_File")
                    .files.values_list("name", flat=True)
                ),
                [case["script"]],
            )
            self.assertEqual(
                list(
                    run.inputs.get(name="ePolyScat_Input_Data")
                    .files.values_list("name", flat=True)
                ),
                [case["runtime_data"]],
            )

    def test_steady_tutorial_sync_is_read_only_and_preserves_timestamps(self):
        tutorial_view = tutorials.ensure_tutorial_view()
        original_view_updated = tutorial_view.updated
        original_run_updates = dict(
            tutorial_view.runs.values_list("pk", "updated")
        )

        with CaptureQueriesContext(connection) as queries:
            synchronized_view = tutorials.ensure_tutorial_view()

        write_sql = [
            query["sql"]
            for query in queries.captured_queries
            if query["sql"].lstrip().upper().startswith(
                ("INSERT", "UPDATE", "DELETE")
            )
        ]
        synchronized_view.refresh_from_db()
        self.assertEqual(synchronized_view.updated, original_view_updated)
        self.assertEqual(
            dict(synchronized_view.runs.values_list("pk", "updated")),
            original_run_updates,
        )
        self.assertEqual(write_sql, [])

    def test_tutorial_database_keys_prevent_duplicate_shared_records(self):
        tutorial_view = tutorials.ensure_tutorial_view()
        tutorial_run = tutorial_view.runs.order_by("pk").first()

        with self.assertRaises(IntegrityError), transaction.atomic():
            models.View.objects.create(
                name="Duplicate Tutorials",
                owner=None,
                type="tutorial",
                tutorial_key=tutorial_view.tutorial_key,
            )
        with self.assertRaises(IntegrityError), transaction.atomic():
            models.Run.objects.create(
                name="Duplicate tutorial run",
                owner=None,
                tutorial_key=tutorial_run.tutorial_key,
            )

    def test_tutorial_view_action_is_visible_with_nested_completed_runs(self):
        for case_path in reversed(tutorials.discover_cases()):
            models.Run.objects.create(
                owner=None,
                name=f"seed-{case_path.name}",
                number=case_path.name,
                filepath=case_path.relative_to(tutorials.PACKAGE_DIR).as_posix(),
            )

        request = self.authenticated_request(
            "get",
            "/epolyscat_django_app/api/views/tutorials/",
        )
        response = views.ViewsViewSet.as_view({"get": "tutorials"})(request)

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.data["name"], "Tutorials")
        self.assertEqual(response.data["type"], "tutorial")
        self.assertTrue(response.data["is_tutorial"])
        self.assertEqual(len(response.data["runs"]), 5)
        self.assertEqual(
            [run_data["name"] for run_data in response.data["runs"]],
            [
                "Test 01",
                "Test 04",
                "Test 08",
                "Test 12",
                "Test 36",
            ],
        )
        self.assertTrue(
            all(run_data["status"] == "COMPLETED" for run_data in response.data["runs"])
        )
        self.assertTrue(
            all(run_data["job_status"] == "COMPLETED" for run_data in response.data["runs"])
        )

    def test_legacy_tutorial_runs_action_synchronizes_the_same_five_runs(self):
        request = self.authenticated_request(
            "get",
            "/epolyscat_django_app/api/runs/tutorial_runs/",
        )
        response = views.RunViewSet.as_view({"get": "tutorial_runs"})(request)

        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            [run_data["name"] for run_data in response.data],
            [
                "Test 01",
                "Test 04",
                "Test 08",
                "Test 12",
                "Test 36",
            ],
        )
        self.assertTrue(all(run_data["is_tutorial"] for run_data in response.data))

    def test_shared_tutorial_does_not_broaden_the_normal_views_list(self):
        tutorials.ensure_tutorial_view()
        ownerless_non_tutorial = models.View.objects.create(
            owner=None,
            name="Unrelated shared view",
            type="default",
        )
        unrelated_keyed_tutorial = models.View.objects.create(
            owner=None,
            name="Unrelated keyed tutorial",
            type="tutorial",
            tutorial_key="unrelated:tutorial",
        )
        user_view = models.View.objects.create(
            owner=self.user,
            name="My view",
            type="user-defined",
        )
        request = self.authenticated_request(
            "get",
            "/epolyscat_django_app/api/views/",
        )

        response = views.ViewsViewSet.as_view({"get": "list"})(request)

        self.assertEqual(response.status_code, 200)
        listed_ids = {view_data["id"] for view_data in response.data["results"]}
        self.assertEqual(listed_ids, {user_view.pk})
        self.assertNotIn(ownerless_non_tutorial.pk, listed_ids)

        retrieve_request = self.authenticated_request(
            "get",
            f"/epolyscat_django_app/api/views/{ownerless_non_tutorial.pk}/",
        )
        retrieve_response = views.ViewsViewSet.as_view({"get": "retrieve"})(
            retrieve_request,
            pk=ownerless_non_tutorial.pk,
        )
        self.assertEqual(retrieve_response.status_code, 404)
        keyed_retrieve_request = self.authenticated_request(
            "get",
            f"/epolyscat_django_app/api/views/{unrelated_keyed_tutorial.pk}/",
        )
        keyed_retrieve_response = views.ViewsViewSet.as_view(
            {"get": "retrieve"}
        )(
            keyed_retrieve_request,
            pk=unrelated_keyed_tutorial.pk,
        )
        self.assertEqual(keyed_retrieve_response.status_code, 404)

    def test_tutorial_status_serialization_never_calls_airavata(self):
        run = tutorials.ensure_tutorial_view().runs.order_by("name").first()
        request = SimpleNamespace(
            user=self.user,
            airavata_client=SimpleNamespace(
                getExperimentStatus=mock.Mock(),
                getJobStatuses=mock.Mock(),
                getJobDetails=mock.Mock(),
            ),
        )

        data = serializers.RunSerializer(run, context={"request": request}).data

        self.assertEqual(data["status"], "COMPLETED")
        self.assertEqual(data["job_status"], "COMPLETED")
        request.airavata_client.getExperimentStatus.assert_not_called()
        request.airavata_client.getJobStatuses.assert_not_called()
        request.airavata_client.getJobDetails.assert_not_called()

    def test_tutorial_view_and_runs_are_read_only_to_authenticated_users(self):
        tutorial_view = tutorials.ensure_tutorial_view()
        run = tutorial_view.runs.first()

        run_request = self.authenticated_request(
            "patch",
            f"/epolyscat_django_app/api/runs/{run.pk}/",
            {"name": "changed"},
        )
        run_response = views.RunViewSet.as_view({"patch": "partial_update"})(
            run_request,
            pk=run.pk,
        )
        view_request = self.authenticated_request(
            "patch",
            f"/epolyscat_django_app/api/views/{tutorial_view.pk}/",
            {"name": "changed", "runIds": [], "overide": True},
        )
        view_response = views.ViewsViewSet.as_view({"patch": "partial_update"})(
            view_request,
            pk=tutorial_view.pk,
        )

        self.assertEqual(run_response.status_code, 403)
        self.assertEqual(view_response.status_code, 403)
        run.refresh_from_db()
        tutorial_view.refresh_from_db()
        self.assertNotEqual(run.name, "changed")
        self.assertEqual(tutorial_view.name, "Tutorials")

    @mock.patch("epolyscat_django_app.views._create_run_user_dir")
    @override_settings(GATEWAY_ID="test-gateway")
    def test_create_cannot_attach_a_user_run_to_the_shared_tutorial_view(
        self,
        _create_run_user_dir,
    ):
        tutorial_view = tutorials.ensure_tutorial_view()
        user_view = models.View.objects.create(
            owner=self.user,
            name="Editable user view",
            type="user-defined",
        )
        create_request = self.authenticated_request(
            "post",
            "/epolyscat_django_app/api/runs/",
            {
                "name": "Attempted tutorial injection",
                "description": "",
                "run_mode": "module",
                "module_application": "ePolyScat",
                "workflow_metadata": {},
                "viewIds": [-1, tutorial_view.pk, user_view.pk],
                "inputs_data": [],
            },
        )
        create_request.authz_token = "token"
        create_request.airavata_client = SimpleNamespace(
            getUserProjects=mock.Mock(
                return_value=[SimpleNamespace(projectID="EPOLYSCAT_app_project")],
            ),
        )

        response = views.RunViewSet.as_view({"post": "create"})(create_request)

        self.assertEqual(response.status_code, 200)
        created_run = models.Run.objects.get(pk=response.data["id"])
        self.assertEqual(created_run.owner, self.user)
        self.assertFalse(created_run.views.filter(pk=tutorial_view.pk).exists())
        self.assertTrue(created_run.views.filter(pk=user_view.pk).exists())
        self.assertFalse(created_run.is_tutorial)

    def test_tutorial_file_catalogs_separate_input_output_and_reference(self):
        run = tutorials.ensure_tutorial_view().runs.get(name="Test 36")
        request = self.authenticated_request(
            "get",
            f"/epolyscat_django_app/api/runs/{run.pk}/input-files/",
        )
        viewset = views.RunViewSet()
        viewset.request = request
        viewset.get_object = mock.Mock(return_value=run)

        input_response = viewset.input_files(request, pk=run.pk)
        output_response = viewset.get_output_files(request, pk=run.pk)

        self.assertEqual(
            [file_data["name"] for file_data in input_response.data],
            ["test36.inp", "test36.gms.dat"],
        )
        self.assertEqual(
            [file_data["name"] for file_data in output_response.data],
            ["test36.ostnd"],
        )
        self.assertTrue(all(file_data["viewable"] for file_data in output_response.data))
        self.assertEqual(tutorials.reference_file_names(run), ["test36.gms.inp"])

        presentation = serializers.RunSerializer(
            run,
            context={"request": request},
        ).data["presentation"]
        self.assertEqual(
            presentation["file_groups"],
            [
                {"label": "Inputs", "files": ["test36.inp", "test36.gms.dat"]},
                {"label": "Output", "files": ["test36.ostnd"]},
                {"label": "Reference", "files": ["test36.gms.inp"]},
            ],
        )

    def test_missing_packaged_artifact_is_reported_as_not_found(self):
        run = tutorials.ensure_tutorial_view().runs.get(name="Test 01")
        preview_request = self.authenticated_request(
            "get",
            f"/epolyscat_django_app/api/runs/{run.pk}/viewables/test01.inp/",
        )
        input_request = self.authenticated_request(
            "get",
            f"/epolyscat_django_app/api/runs/{run.pk}/input-files/",
        )
        output_request = self.authenticated_request(
            "get",
            f"/epolyscat_django_app/api/runs/{run.pk}/get_output_files/",
        )

        with mock.patch.object(
            tutorials,
            "load_catalog",
            side_effect=ValueError("declared artifact is missing"),
        ):
            preview_response = views.RunViewSet.as_view(
                {"get": "show_viewable"}
            )(
                preview_request,
                pk=run.pk,
                filename="test01.inp",
            )
            input_response = views.RunViewSet.as_view({"get": "input_files"})(
                input_request,
                pk=run.pk,
            )
            output_response = views.RunViewSet.as_view(
                {"get": "get_output_files"}
            )(
                output_request,
                pk=run.pk,
            )

        self.assertEqual(preview_response.status_code, 404)
        self.assertEqual(input_response.status_code, 404)
        self.assertEqual(output_response.status_code, 404)

    def test_tutorial_text_preview_reads_only_packaged_case_files(self):
        run = tutorials.ensure_tutorial_view().runs.get(name="Test 36")
        request = self.authenticated_request(
            "get",
            f"/epolyscat_django_app/api/runs/{run.pk}/viewables/test36.inp/",
        )
        viewset = views.RunViewSet()
        viewset.request = request
        viewset.get_object = mock.Mock(return_value=run)

        response = viewset.show_viewable(request, pk=run.pk, filename="test36.inp")
        contents = b"".join(response.streaming_content).decode("utf-8")

        self.assertIn("Convert", contents)

        output_response = viewset.show_viewable(
            request,
            pk=run.pk,
            filename="test36.ostnd",
        )
        output_contents = b"".join(output_response.streaming_content).decode(
            "utf-8",
            errors="replace",
        )
        self.assertTrue(output_contents.strip())

        reference_response = viewset.show_viewable(
            request,
            pk=run.pk,
            filename="test36.gms.inp",
        )
        reference_contents = b"".join(reference_response.streaming_content).decode(
            "utf-8",
            errors="replace",
        )
        self.assertTrue(reference_contents.strip())

        with self.assertRaises(FileNotFoundError):
            tutorials.open_case_file(run, "../test36.inp")

    @mock.patch(
        "epolyscat_django_app.views._save_run_user_file",
        return_value=SimpleNamespace(productUri="clone-product-uri"),
    )
    @mock.patch("epolyscat_django_app.views._create_run_user_dir")
    @override_settings(GATEWAY_ID="test-gateway")
    def test_hydrated_tutorial_clone_is_user_owned_and_editable(
        self,
        _create_run_user_dir,
        _save_run_user_file,
    ):
        tutorial_run = tutorials.ensure_tutorial_view().runs.get(
            name="Test 36",
        )
        with tutorials.open_case_file(tutorial_run, "test36.inp") as input_stream:
            script_contents = input_stream.read().decode("utf-8")
        with tutorials.open_case_file(tutorial_run, "test36.gms.dat") as data_stream:
            runtime_data_contents = data_stream.read().decode("utf-8")

        create_request = self.authenticated_request(
            "post",
            "/epolyscat_django_app/api/runs/",
            {
                "name": f"Clone of {tutorial_run.name}",
                "description": tutorial_run.description,
                "run_mode": tutorial_run.run_mode,
                "module_application": tutorial_run.module_application,
                "workflow_metadata": {},
                "viewIds": [],
                "inputs_data": [
                    {
                        "type": "radio input",
                        "name": "Calculation_Type",
                        "value": "MODULE",
                    },
                    {
                        "type": "radio input",
                        "name": "EPOLYSCAT_Application_Module",
                        "value": "ePolyScat",
                    },
                    {
                        "type": "files",
                                "name": "ePolyscat_Input_File",
                        "files": [
                            {
                                "name": "test36.inp",
                                "contents": script_contents,
                                "isPlaintext": True,
                                "deleted": False,
                            },
                        ],
                    },
                    {
                        "type": "files",
                        "name": "ePolyScat_Input_Data",
                        "files": [
                            {
                                "name": "test36.gms.dat",
                                "contents": runtime_data_contents,
                                "isPlaintext": True,
                                "deleted": False,
                            },
                        ],
                    },
                ],
            },
        )
        create_request.authz_token = "token"
        create_request.airavata_client = SimpleNamespace(
            getUserProjects=mock.Mock(
                return_value=[SimpleNamespace(projectID="EPOLYSCAT_app_project")],
            ),
        )

        create_response = views.RunViewSet.as_view({"post": "create"})(
            create_request,
        )

        self.assertEqual(create_response.status_code, 200)
        clone = models.Run.objects.get(pk=create_response.data["id"])
        self.assertEqual(clone.owner, self.user)
        self.assertFalse(clone.is_tutorial)
        self.assertFalse(clone.views.filter(type="tutorial").exists())
        self.assertEqual(
            set(clone.inputs.filter(type="files").values_list("name", flat=True)),
            {"ePolyscat_Input_File", "ePolyScat_Input_Data"},
        )
        self.assertEqual(
            clone.inputs.get(name="ePolyscat_Input_File")
            .files.get(name="test36.inp")
            .data_product_uri,
            "clone-product-uri",
        )
        self.assertEqual(
            clone.inputs.get(name="ePolyScat_Input_Data")
            .files.get(name="test36.gms.dat")
            .data_product_uri,
            "clone-product-uri",
        )
        self.assertEqual(_save_run_user_file.call_count, 2)

        update_request = self.authenticated_request(
            "patch",
            f"/epolyscat_django_app/api/runs/{clone.pk}/",
            {"name": "My editable tutorial clone", "inputs_data": []},
        )
        update_response = views.RunViewSet.as_view({"patch": "partial_update"})(
            update_request,
            pk=clone.pk,
        )

        self.assertEqual(update_response.status_code, 200)
        clone.refresh_from_db()
        self.assertEqual(clone.name, "My editable tutorial clone")

    def test_current_mutating_apis_reject_tutorial_sources_before_external_calls(self):
        tutorial_view = tutorials.ensure_tutorial_view()
        run = tutorial_view.runs.first()
        snapshot = (run.name, run.deleted, run.updated, tuple(run.views.values_list('pk', flat=True)))
        actions = [
            ('put', 'update', views.RunViewSet, run.pk),
            ('delete', 'destroy', views.RunViewSet, run.pk),
            ('post', 'submit', views.RunViewSet, run.pk),
            ('post', 'resubmit', views.RunViewSet, run.pk),
            ('put', 'cancel', views.RunViewSet, run.pk),
            ('patch', 'change_notification_settings', views.RunViewSet, run.pk),
            ('post', 'workflow_continuation', views.RunViewSet, run.pk),
            ('delete', 'destroy', views.ViewsViewSet, tutorial_view.pk),
        ]
        for method, action, viewset, pk in actions:
            with self.subTest(action=action):
                request = self.authenticated_request(method, '/isolated-tutorial-check/', {})
                request.airavata_client = mock.Mock()
                response = viewset.as_view({method: action})(request, pk=pk)
                self.assertEqual(response.status_code, 403)
                self.assertEqual(request.airavata_client.mock_calls, [])
        run.refresh_from_db()
        self.assertEqual(snapshot, (run.name, run.deleted, run.updated, tuple(run.views.values_list('pk', flat=True))))

    def test_key_migration_preserves_legacy_rows_and_repeat_sync_is_separate(self):
        import importlib
        from django.apps import apps
        migration = importlib.import_module('epolyscat_django_app.migrations.0028_tutorial_keys')
        first_view = models.View.objects.create(name='Legacy Tutorials', type='tutorial', owner=None)
        second_view = models.View.objects.create(name='Legacy duplicate', type='tutorial', owner=None)
        first_run = models.Run.objects.create(name='Legacy Test 01', owner=None, filepath='data/epolyscat/manual_tutorials/test01')
        second_run = models.Run.objects.create(name='Legacy duplicate', owner=None, filepath=first_run.filepath)
        first_view.runs.add(first_run)
        second_view.runs.add(second_run)
        user_run = models.Run.objects.create(name='Existing user run', owner=self.user, filepath=first_run.filepath)
        before = list(models.Run.objects.values())
        before_views = list(models.View.objects.values())
        for _ in range(2):
            migration.assign_tutorial_keys(apps, None)
        migration.clear_tutorial_keys(apps, None)
        self.assertEqual(list(models.Run.objects.values()), before)
        self.assertEqual(list(models.View.objects.values()), before_views)
        canonical = tutorials.ensure_tutorial_view()
        self.assertNotIn(canonical.pk, [first_view.pk, second_view.pk])
        self.assertEqual(canonical.runs.count(), 5)
        self.assertEqual(tutorials.ensure_tutorial_view().pk, canonical.pk)
        self.assertEqual(list(models.Run.objects.filter(pk__in=[first_run.pk,second_run.pk,user_run.pk]).values()), before)
        self.assertTrue(first_view.runs.filter(pk=first_run.pk).exists())
        self.assertTrue(second_view.runs.filter(pk=second_run.pk).exists())

    def test_legacy_personal_relations_and_history_are_not_adopted(self):
        admin = get_user_model().objects.create_user(username='synthetic-legacy-admin')
        root = models.RunsRoot.objects.create(root='Legacy root',owner=admin)
        experiment = models.Experiment.objects.create(name='Legacy experiment',description='',owner=admin)
        view = models.View.objects.create(name='Legacy Tutorials',owner=None,type='tutorial',tutorial_key='epolyscat:manual-tutorials')
        run = models.Run.objects.create(name='Legacy Test 01',owner=None,root=root,experiment=experiment,filepath='data/epolyscat/manual_tutorials/test01',directory='admin/shared/tutorials',airavata_project_id='legacy-project',inpc_data_product_uri='synthetic-admin-uri',tutorial_key='epolyscat:manual-tutorial:test01',workflow_metadata={'tutorialCaseId':'test01','tutorialStandalone':True})
        view.runs.add(run)
        execution = models.RemoteExecution.objects.create(run=run,airavata_experiment_id='legacy-experiment',airavata_experiment_status='COMPLETED')
        source_input = models.Input.objects.create(run=run,type='files',name='ePolyscat_Input_File')
        source_file = models.File.objects.create(input=source_input,name='test01.inp',data_product_uri='synthetic-admin-uri')
        before = models.Run.objects.filter(pk=run.pk).values().get()
        view_before = models.View.objects.filter(pk=view.pk).values().get()
        canonical = tutorials.ensure_tutorial_view()
        self.assertNotEqual(canonical.pk,view.pk)
        self.assertEqual(models.Run.objects.filter(pk=run.pk).values().get(),before)
        self.assertEqual(models.View.objects.filter(pk=view.pk).values().get(),view_before)
        self.assertTrue(view.runs.filter(pk=run.pk).exists())
        self.assertTrue(models.RemoteExecution.objects.filter(pk=execution.pk).exists())
        source_file.refresh_from_db()
        self.assertEqual(source_file.data_product_uri,'synthetic-admin-uri')
        request = self.authenticated_request('get',f'/runs/{run.pk}/viewables/test01.inp/')
        with mock.patch('epolyscat_django_app.views.user_storage.open_file') as open_remote:
            response = views.RunViewSet.as_view({'get':'show_viewable'})(request,pk=run.pk,filename='test01.inp')
        self.assertEqual(response.status_code,404)
        open_remote.assert_not_called()
        admin.delete()
        canonical.refresh_from_db()
        self.assertEqual(canonical.runs.count(),5)
        self.assertEqual(tutorials.ensure_tutorial_view().pk,canonical.pk)
        response = views.ViewsViewSet.as_view({'get':'tutorials'})(self.authenticated_request('get','/views/tutorials/'))
        self.assertEqual(response.status_code,200)
        self.assertEqual(len(response.data['runs']),5)
        for source in canonical.runs.all():
            self.assertIsNone(source.owner_id)
            self.assertIsNone(source.root_id)
            self.assertIsNone(source.experiment_id)
            self.assertIsNone(source.airavata_project_id)
            self.assertEqual(source.directory,'')
            self.assertFalse(source.executions.exists())
        source=canonical.runs.get(name='Test 01')
        preview = views.RunViewSet.as_view({'get':'show_viewable'})(self.authenticated_request('get','/preview/'),pk=source.pk,filename='test01.inp')
        self.assertEqual(preview.status_code,200)
        body=b''.join(preview.streaming_content).decode('utf-8')
        self.assertIn('Convert',body)
        # This existing end-to-end API clone regression creates a new user run
        # from package bytes, stores both inputs, and verifies editable ownership.
        self.test_hydrated_tutorial_clone_is_user_owned_and_editable()

    def test_database_rejects_personal_links_on_packaged_sources(self):
        canonical=tutorials.ensure_tutorial_view()
        run=canonical.runs.first()
        root=models.RunsRoot.objects.create(root='User root',owner=self.user)
        experiment=models.Experiment.objects.create(name='User experiment',description='',owner=self.user)
        for attrs in [{'owner':self.user},{'root':root},{'experiment':experiment},{'directory':'personal/path'},{'airavata_project_id':'personal-project'},{'inpc_data_product_uri':'personal-uri'}]:
            with self.subTest(attrs=attrs),self.assertRaises(IntegrityError),transaction.atomic():
                models.Run.objects.filter(pk=run.pk).update(**attrs)
        with self.assertRaises(IntegrityError),transaction.atomic():
            models.View.objects.filter(pk=canonical.pk).update(owner=self.user)

    def test_missing_catalog_fails_explicitly_without_remote_fallback(self):
        tutorials.ensure_tutorial_view()
        with mock.patch.object(tutorials,'load_catalog',side_effect=ValueError('missing packaged artifact')),mock.patch('epolyscat_django_app.views.user_storage.open_file') as open_remote:
            response=views.ViewsViewSet.as_view({'get':'tutorials'})(self.authenticated_request('get','/views/tutorials/'))
        self.assertEqual(response.status_code,503)
        self.assertIn('no personal storage fallback',str(response.data['detail']))
        open_remote.assert_not_called()

    def test_run_listing_uses_canonical_view_even_when_legacy_view_is_older(self):
        legacy=models.View.objects.create(name='Older legacy Tutorials',type='tutorial',owner=None)
        canonical=tutorials.ensure_tutorial_view()
        self.assertLess(legacy.pk,canonical.pk)
        response=views.RunViewSet.as_view({'get':'list'})(self.authenticated_request('get',f'/runs/?viewId={canonical.pk}'))
        self.assertEqual(response.status_code,200)
        self.assertEqual(len(response.data['results']),5)
        self.assertTrue(all(row['is_tutorial'] and row['status']=='COMPLETED' for row in response.data['results']))

    def test_packaged_source_rejects_new_execution_without_removing_legacy_history(self):
        from django.core.exceptions import ValidationError
        canonical=tutorials.ensure_tutorial_view().runs.first()
        with self.assertRaises(ValidationError):
            models.RemoteExecution.objects.create(run=canonical,airavata_experiment_id='forbidden-canonical-execution')
        self.assertFalse(canonical.executions.exists())
        legacy=models.Run.objects.create(name='Legacy execution source',owner=None,tutorial_key='epolyscat:manual-tutorial:test01')
        execution=models.RemoteExecution.objects.create(run=legacy,airavata_experiment_id='retained-legacy-execution')
        tutorials.ensure_tutorial_view()
        self.assertTrue(models.RemoteExecution.objects.filter(pk=execution.pk).exists())
