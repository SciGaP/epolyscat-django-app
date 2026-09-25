from pathlib import Path

from django.test import SimpleTestCase


APP_DIR = Path(__file__).resolve().parents[1]
THRIFT_RUNTIME_MODULES = (
    "models.py",
    "serializers.py",
    "scientific_output_service.py",
    "views.py",
)


class ThriftRuntimeContractTests(SimpleTestCase):
    def test_runtime_uses_only_thrift_airavata_dependencies(self):
        runtime_source = "\n".join(
            (APP_DIR / filename).read_text()
            for filename in THRIFT_RUNTIME_MODULES
        )

        self.assertNotIn("airavata_grpc", runtime_source)
        self.assertNotIn("airavata_sdk.generated", runtime_source)
        self.assertNotIn("airavata_channel", runtime_source)
        self.assertNotIn("request.airavata.", runtime_source)
        self.assertIn("request.airavata_client", runtime_source)
        self.assertIn("request.authz_token", runtime_source)

    def test_grpc_compatibility_module_is_removed(self):
        self.assertFalse((APP_DIR / "airavata_grpc.py").exists())
