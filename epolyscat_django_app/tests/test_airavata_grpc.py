from enum import IntEnum
from unittest import TestCase

from epolyscat_django_app.airavata_grpc import _thrift_enum_compat


class ThriftEnumCompatibilityTests(TestCase):
    def test_supports_python_int_enum_generated_by_newer_sdk(self):
        class ModernExperimentState(IntEnum):
            CREATED = 0
            COMPLETED = 7
            FAILED = 8

        compat = _thrift_enum_compat(ModernExperimentState)

        self.assertEqual(0, compat.CREATED)
        self.assertEqual(7, compat.COMPLETED)
        self.assertEqual("COMPLETED", compat(compat.COMPLETED).name)
        self.assertEqual(7, compat["COMPLETED"].value)
        self.assertEqual("FAILED", compat.Name(8))
        self.assertEqual(0, compat.Value("CREATED"))

    def test_preserves_legacy_thrift_enum_maps(self):
        class LegacyExperimentState:
            _VALUES_TO_NAMES = {0: "CREATED", 7: "COMPLETED"}
            _NAMES_TO_VALUES = {"CREATED": 0, "COMPLETED": 7}

        compat = _thrift_enum_compat(LegacyExperimentState)

        self.assertEqual(0, compat.CREATED)
        self.assertEqual("COMPLETED", compat(7).name)
        self.assertEqual(7, compat.Value("COMPLETED"))
