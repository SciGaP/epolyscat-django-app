from enum import IntEnum
from unittest import TestCase

from epolyscat_django_app.airavata_grpc import (
    _thrift_enum_compat,
    replace_model_list,
)


class _RepeatedValues:
    def __init__(self, values):
        self.values = list(values)

    def __delitem__(self, key):
        del self.values[key]

    def extend(self, values):
        self.values.extend(values)


class _MessageWithRepeatedValues:
    def __init__(self, values):
        self.email_addresses = _RepeatedValues(values)


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


class ProtobufFieldHelpersTests(TestCase):
    def test_replace_model_list_clears_existing_repeated_values(self):
        message = _MessageWithRepeatedValues(["old@example.org"])

        replace_model_list(
            message,
            "email_addresses",
            ["new@example.org"],
        )

        self.assertEqual(
            message.email_addresses.values,
            ["new@example.org"],
        )
