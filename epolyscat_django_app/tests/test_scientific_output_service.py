from unittest import mock

from django.contrib.auth import get_user_model
from django.core.cache import cache
from django.test import RequestFactory, TestCase

from epolyscat_django_app import models, scientific_output_service


class ScientificOutputCacheTests(TestCase):
    def setUp(self):
        cache.clear()
        self.user = get_user_model().objects.create_user(
            username="scientific-output-cache"
        )
        experiment = models.Experiment.objects.create(
            name="experiment",
            owner=self.user,
        )
        root = models.RunsRoot.objects.create(
            root="root",
            owner=self.user,
        )
        self.run = models.Run.objects.create(
            name="run",
            owner=self.user,
            root=root,
            number="0001",
            experiment=experiment,
            module_application="Gaussian16",
        )
        models.RemoteExecution.objects.create(
            run=self.run,
            airavata_experiment_id="experiment-id",
            airavata_experiment_status="COMPLETED",
        )
        self.request = RequestFactory().get("/")
        self.request.user = self.user

    def tearDown(self):
        cache.clear()

    @mock.patch(
        "epolyscat_django_app.scientific_output_service."
        "user_storage.list_experiment_dir"
    )
    def test_completed_output_listing_is_reused(self, list_experiment_dir):
        list_experiment_dir.return_value = (
            [],
            [{"name": "gaussian.log"}],
        )

        first = scientific_output_service.output_files_for_run(
            self.request,
            self.run,
        )
        second = scientific_output_service.output_files_for_run(
            self.request,
            self.run,
        )

        self.assertEqual(first, second)
        self.assertEqual(list_experiment_dir.call_count, 2)

    @mock.patch(
        "epolyscat_django_app.scientific_output_service.cache.set"
    )
    @mock.patch(
        "epolyscat_django_app.scientific_output_service."
        "user_storage.list_experiment_dir"
    )
    def test_running_output_listing_uses_short_cache_ttl(
        self,
        list_experiment_dir,
        cache_set,
    ):
        execution = self.run.latest_execution
        execution.airavata_experiment_status = "EXECUTING"
        execution.save(update_fields=["airavata_experiment_status"])
        list_experiment_dir.return_value = ([], [])

        scientific_output_service.output_files_for_run(
            self.request,
            self.run,
        )

        self.assertEqual(
            cache_set.call_args.kwargs["timeout"],
            scientific_output_service.RUNNING_OUTPUT_CACHE_SECONDS,
        )
