import importlib
import importlib.util

from django.contrib.auth import get_user_model
from django.db import models as django_models
from django.test import RequestFactory
from django.test import SimpleTestCase
from django.test import TestCase

from epolyscat_django_app import models, serializers
from epolyscat_django_app.workflow_continuation import (
    classify_run,
    continuation_eligibility,
)


class WorkflowContinuationModelTests(SimpleTestCase):
    def test_run_references_imported_workflow_source_without_reparenting_it(self):
        field = models.Run._meta.get_field("workflow_source_run")

        self.assertIsInstance(field, django_models.ForeignKey)
        self.assertIs(field.remote_field.model, models.Run)
        self.assertIs(field.remote_field.on_delete, django_models.PROTECT)
        self.assertEqual(field.remote_field.related_name, "workflow_continuations")
        self.assertTrue(field.null)
        self.assertTrue(field.blank)

    def test_run_serializer_exposes_workflow_source_run(self):
        self.assertIn("workflow_source_run", serializers.RunSerializer.Meta.fields)


class WorkflowContinuationDomainTests(SimpleTestCase):
    def test_workflow_continuation_domain_exposes_run_classifier(self):
        module_name = "epolyscat_django_app.workflow_continuation"

        self.assertIsNotNone(importlib.util.find_spec(module_name))
        module = importlib.import_module(module_name)
        self.assertTrue(callable(module.classify_run))

    def test_workflow_continuation_domain_exposes_eligibility_check(self):
        module = importlib.import_module(
            "epolyscat_django_app.workflow_continuation"
        )

        self.assertTrue(hasattr(module, "continuation_eligibility"))
        self.assertTrue(callable(module.continuation_eligibility))

    def test_classifies_data_generation_modules(self):
        for application in ("Gaussian16", "OpenMolcas"):
            with self.subTest(application=application):
                run = models.Run(run_mode="module", module_application=application)

                self.assertEqual(
                    classify_run(run),
                    {
                        "source_stage": "Data_Gen",
                        "source_application": application,
                        "next_stage": "ePolyScat_Run",
                    },
                )

    def test_classifies_epolyscat_module(self):
        run = models.Run(run_mode="module", module_application="ePolyScat")

        self.assertEqual(
            classify_run(run),
            {
                "source_stage": "ePolyScat_Run",
                "source_application": "ePolyScat",
                "next_stage": "Analysis",
            },
        )

    def test_classifies_analysis_utility_as_terminal(self):
        run = models.Run(run_mode="utility", utility_application="CnvLinFull")

        self.assertEqual(
            classify_run(run),
            {
                "source_stage": "Analysis",
                "source_application": "CnvLinFull",
                "next_stage": None,
            },
        )

    def test_classifies_workflow_child_from_its_stage(self):
        run = models.Run(
            run_mode="workflow",
            workflow_stage="Data_Gen",
            workflow_application="OpenMolcas",
        )

        self.assertEqual(
            classify_run(run),
            {
                "source_stage": "Data_Gen",
                "source_application": "OpenMolcas",
                "next_stage": "ePolyScat_Run",
            },
        )

    def test_rejects_workflow_parent_and_unknown_module(self):
        workflow_parent = models.Run(
            run_mode="workflow",
            workflow_metadata={"isWorkflowPlan": True},
        )
        unknown_module = models.Run(
            run_mode="module",
            module_application="UnknownApplication",
        )

        self.assertIsNone(classify_run(workflow_parent))
        self.assertIsNone(classify_run(unknown_module))


class WorkflowContinuationLegacyRunTests(TestCase):
    def test_classifies_legacy_run_selector_inputs_when_new_fields_are_blank(self):
        run = models.Run.objects.create(name="legacy utility")
        models.Input.objects.create(
            run=run,
            type="radio input",
            name="Calculation_Type",
            value="UTILITY",
        )
        models.Input.objects.create(
            run=run,
            type="radio input",
            name="Application_Utility",
            value="NRFPAD",
        )

        self.assertEqual(
            classify_run(run),
            {
                "source_stage": "Analysis",
                "source_application": "NRFPAD",
                "next_stage": None,
            },
        )


class WorkflowContinuationEligibilityTests(TestCase):
    def setUp(self):
        self.user = get_user_model().objects.create_user(username="workflow-user")
        self.request = RequestFactory().get("/api/runs/1/workflow_continuation/")
        self.request.user = self.user

    def create_run(self, application="ePolyScat", status=None, **overrides):
        values = {
            "name": "source run",
            "owner": self.user,
            "run_mode": "module",
            "module_application": application,
        }
        values.update(overrides)
        run = models.Run.objects.create(**values)
        if status:
            models.RemoteExecution.objects.create(
                run=run,
                airavata_experiment_id=f"experiment-{run.id}",
                airavata_experiment_status=status,
            )
        return run

    def test_completed_supported_run_is_eligible(self):
        run = self.create_run(status="COMPLETED")

        result = continuation_eligibility(run, self.request)

        self.assertIn("eligible", result)
        self.assertTrue(result["eligible"])
        self.assertEqual(result["reason"], "")
        self.assertEqual(result["source_stage"], "ePolyScat_Run")
        self.assertEqual(result["source_application"], "ePolyScat")
        self.assertEqual(result["next_stage"], "Analysis")

    def test_run_without_execution_is_not_eligible(self):
        run = self.create_run()

        result = continuation_eligibility(run, self.request)

        self.assertIn("eligible", result)
        self.assertFalse(result["eligible"])
        self.assertEqual(result["reason"], "no_execution")

    def test_failed_run_is_not_eligible(self):
        run = self.create_run(status="FAILED")

        result = continuation_eligibility(run, self.request)

        self.assertIn("eligible", result)
        self.assertFalse(result["eligible"])
        self.assertEqual(result["reason"], "run_not_completed")

    def test_unknown_run_type_is_not_eligible(self):
        run = self.create_run(application="Unknown", status="COMPLETED")

        result = continuation_eligibility(run, self.request)

        self.assertIn("eligible", result)
        self.assertFalse(result["eligible"])
        self.assertEqual(result["reason"], "unsupported_run_type")

    def test_analysis_utility_is_terminal(self):
        run = self.create_run(
            application="",
            status="COMPLETED",
            run_mode="utility",
            utility_application="CnvMath",
        )

        result = continuation_eligibility(run, self.request)

        self.assertIn("eligible", result)
        self.assertFalse(result["eligible"])
        self.assertEqual(result["reason"], "terminal_stage")
