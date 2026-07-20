"""Compatibility layer for Thrift and gRPC Airavata portal runtimes.

The dev portal currently runs the gRPC Airavata SDK, while local and older
portal installs may still have the Thrift SDK. Keep the rest of this app using
the old small surface: enum objects with ``.name``/``.value``, ``user_storage``,
``experiment_util``, and ``remoteapi``.
"""

from types import SimpleNamespace


_LEGACY_FIELD_NAMES = {
    "application_interface_id": "applicationInterfaceId",
    "application_inputs": "applicationInputs",
    "application_modules": "applicationModules",
    "application_outputs": "applicationOutputs",
    "app_module_id": "appModuleId",
    "compute_host_id": "computeHostId",
    "computational_resource_scheduling": "computationalResourceScheduling",
    "email_addresses": "emailAddresses",
    "execution_id": "executionId",
    "executable_path": "executablePath",
    "experiment_id": "experimentId",
    "experiment_inputs": "experimentInputs",
    "experiment_name": "experimentName",
    "experiment_outputs": "experimentOutputs",
    "experiment_status": "experimentStatus",
    "gateway_id": "gatewayId",
    "group_resource_profile_id": "groupResourceProfileId",
    "host_name": "hostName",
    "input_order": "inputOrder",
    "is_required": "isRequired",
    "job_id": "jobId",
    "job_state": "jobState",
    "product_uri": "productUri",
    "project_id": "projectId",
    "required_to_added_to_command_line": "requiredToAddedToCommandLine",
    "resource_host_id": "resourceHostId",
    "share_experiment_publicly": "shareExperimentPublicly",
    "time_of_state_change": "timeOfStateChange",
    "total_cpu_count": "totalCPUCount",
    "total_physical_memory": "totalPhysicalMemory",
    "user_configuration_data": "userConfigurationData",
    "user_name": "userName",
    "wall_time_limit": "wallTimeLimit",
}


def _camel_case(name):
    if name in _LEGACY_FIELD_NAMES:
        return _LEGACY_FIELD_NAMES[name]
    head, *tail = name.split("_")
    return head + "".join(part.capitalize() for part in tail)


def model_field(message, name):
    """Read a protobuf field while keeping local legacy-model tests usable."""
    if hasattr(message, name):
        return getattr(message, name)
    return getattr(message, _camel_case(name))


def set_model_field(message, name, value):
    field_name = name if hasattr(message, name) else _camel_case(name)
    setattr(message, field_name, value)


def replace_model_list(message, name, values):
    field_name = name if hasattr(message, name) else _camel_case(name)
    target = getattr(message, field_name, None)
    if hasattr(target, "extend") and not isinstance(target, list):
        del target[:]
        target.extend(values)
    else:
        setattr(message, field_name, list(values))


def copy_model_field(message, name, value):
    field_name = name if hasattr(message, name) else _camel_case(name)
    target = getattr(message, field_name, None)
    if hasattr(target, "CopyFrom"):
        target.CopyFrom(value)
    else:
        setattr(message, field_name, value)


def create_model(model_class, **fields):
    try:
        return model_class(**fields)
    except TypeError:
        return model_class(
            **{_camel_case(name): value for name, value in fields.items()}
        )


class _EnumMember:
    def __init__(self, name, value):
        self.name = name
        self.value = value

    def __int__(self):
        return self.value


def _thrift_enum_compat(enum_cls):
    values_to_names = dict(getattr(enum_cls, "_VALUES_TO_NAMES", {}))
    names_to_values = dict(getattr(enum_cls, "_NAMES_TO_VALUES", {}))

    enum_members = getattr(enum_cls, "__members__", {})
    if enum_members and not names_to_values:
        names_to_values = {
            name: int(member.value) for name, member in enum_members.items()
        }
    if names_to_values and not values_to_names:
        values_to_names = {}
        for name, value in names_to_values.items():
            values_to_names.setdefault(value, name)
    if values_to_names and not names_to_values:
        names_to_values = {
            name: value for value, name in values_to_names.items()
        }

    class CompatEnum:
        _VALUES_TO_NAMES = values_to_names
        _NAMES_TO_VALUES = names_to_values

        def __new__(cls, value):
            return _EnumMember(values_to_names.get(value, str(value)), value)

        def __class_getitem__(cls, name):
            return _EnumMember(name, names_to_values[name])

        @staticmethod
        def Name(value):
            return values_to_names[value]

        @staticmethod
        def Value(name):
            return names_to_values[name]

    for name, value in names_to_values.items():
        setattr(CompatEnum, name, value)
    return CompatEnum


try:
    from airavata_sdk.generated.org.apache.airavata.model.application.io import (
        application_io_pb2 as _io_pb2,
    )
    from airavata_sdk.generated.org.apache.airavata.model.experiment import (
        experiment_pb2 as _exp_pb2,
    )
    from airavata_sdk.generated.org.apache.airavata.model.scheduling import (
        scheduling_pb2 as _sched_pb2,
    )
    from airavata_sdk.generated.org.apache.airavata.model.status import (
        status_pb2 as _status_pb2,
    )
    from airavata_sdk.generated.org.apache.airavata.model.workspace import (
        workspace_pb2 as _ws_pb2,
    )
except ModuleNotFoundError:
    from airavata.model.application.io.ttypes import DataType
    from airavata.model.experiment.ttypes import (
        ExperimentModel,
        UserConfigurationDataModel,
    )
    from airavata.model.scheduling.ttypes import (
        ComputationalResourceSchedulingModel,
    )
    from airavata.model.status.ttypes import ExperimentState as _ExperimentState
    from airavata.model.workspace.ttypes import Project
    from airavata_django_portal_sdk import experiment_util, remoteapi, user_storage

    ExperimentState = _thrift_enum_compat(_ExperimentState)

    def django_user(request):
        from django.contrib.auth import get_user_model

        user = getattr(request, "user", None)
        username = getattr(user, "username", None)
        if not username:
            return user
        User = get_user_model()
        if isinstance(user, User):
            return user
        obj, _ = User.objects.get_or_create(
            username=username,
            defaults={
                "email": getattr(user, "email", "") or "",
                "first_name": getattr(user, "first_name", "") or "",
                "last_name": getattr(user, "last_name", "") or "",
            },
        )
        return obj

else:
    ExperimentModel = _exp_pb2.ExperimentModel
    UserConfigurationDataModel = _exp_pb2.UserConfigurationDataModel
    ComputationalResourceSchedulingModel = (
        _sched_pb2.ComputationalResourceSchedulingModel
    )
    Project = _ws_pb2.Project

    def _enum_compat(enum_wrapper, module, *, strip_prefix=""):
        values_to_names = {}
        names_to_values = {}
        for value in enum_wrapper.DESCRIPTOR.values:
            public_name = value.name
            if strip_prefix and public_name.startswith(strip_prefix):
                public_name = public_name[len(strip_prefix) :]
            values_to_names[value.number] = public_name
            names_to_values[public_name] = value.number

        class CompatEnum:
            DESCRIPTOR = enum_wrapper.DESCRIPTOR

            def __new__(cls, value):
                return _EnumMember(values_to_names.get(value, str(value)), value)

            def __class_getitem__(cls, name):
                return _EnumMember(name, names_to_values[name])

            @staticmethod
            def Name(value):
                return values_to_names[value]

            @staticmethod
            def Value(name):
                return names_to_values[name]

        for name, value in names_to_values.items():
            setattr(CompatEnum, name, value)
        for name in dir(module):
            if name.isupper() and not hasattr(CompatEnum, name):
                value = getattr(module, name)
                if isinstance(value, int):
                    setattr(CompatEnum, name, value)
        return CompatEnum

    DataType = _enum_compat(_io_pb2.DataType, _io_pb2)
    ExperimentState = _enum_compat(
        _status_pb2.ExperimentState,
        _status_pb2,
        strip_prefix="EXPERIMENT_STATE_",
    )

    def _client(request):
        return request.airavata

    def django_user(request):
        from django.contrib.auth import get_user_model

        user = getattr(request, "user", None)
        username = getattr(user, "username", None)
        if not username:
            return user
        User = get_user_model()
        if isinstance(user, User):
            return user
        obj, _ = User.objects.get_or_create(
            username=username,
            defaults={
                "email": getattr(user, "email", "") or "",
                "first_name": getattr(user, "first_name", "") or "",
                "last_name": getattr(user, "last_name", "") or "",
            },
        )
        return obj

    def _path(path=None, *, data_product_uri=None, data_product=None, dir_names=None):
        if path is None and dir_names is not None:
            path = "/".join(dir_names)
        if path is None:
            path = data_product_uri
        if path is None and data_product is not None:
            path = getattr(data_product, "productUri", None) or getattr(
                data_product, "value", None
            )
        return path or ""

    def _data_product(value, name=""):
        if hasattr(value, "productUri"):
            return value
        if hasattr(value, "product_uri"):
            return SimpleNamespace(productUri=value.product_uri, name=name)
        return SimpleNamespace(productUri=value, name=name)

    class _UserStorage:
        @staticmethod
        def exists(request, path=None, **kwargs):
            return _client(request).storage.file_exists(_path(path, **kwargs))

        @staticmethod
        def user_file_exists(request, path=None, **kwargs):
            resolved_path = _path(path, **kwargs)
            try:
                if _client(request).storage.file_exists(resolved_path):
                    return resolved_path
                return None
            except Exception:
                return None

        @staticmethod
        def dir_exists(request, path=None, **kwargs):
            return _client(request).storage.dir_exists(_path(path, **kwargs))

        @staticmethod
        def open_file(request, path=None, **kwargs):
            import io

            return io.BytesIO(
                _client(request).storage.download_file(_path(path, **kwargs))
            )

        @staticmethod
        def save(request, path, file, name=None, **kwargs):
            content = file.read() if hasattr(file, "read") else file
            name = name or getattr(file, "name", None) or path.rsplit("/", 1)[-1]
            upload_path = f"{path.rstrip('/')}/{name}" if path else name
            return _data_product(
                _client(request).storage.upload_file(
                    path=upload_path,
                    content=content,
                    name=name,
                ),
                name=name,
            )

        @staticmethod
        def save_input_file(request, file, name=None, **kwargs):
            content = file.read() if hasattr(file, "read") else file
            name = name or getattr(file, "name", "input")
            return _data_product(
                _client(request).storage.upload_file(
                    path=name,
                    content=content,
                    name=name,
                ),
                name=name,
            )

        @staticmethod
        def create_user_dir(request, path=None, **kwargs):
            return _client(request).storage.create_dir(_path(path, **kwargs))

        @staticmethod
        def delete(request, path=None, **kwargs):
            return _client(request).storage.delete_file(_path(path, **kwargs))

        @staticmethod
        def delete_dir(request, path=None, **kwargs):
            return _client(request).storage.delete_dir(_path(path, **kwargs))

        @staticmethod
        def listdir(request, path=None, **kwargs):
            return _client(request).storage.list_dir(_path(path, **kwargs))

        @staticmethod
        def list_experiment_dir(request, experiment_id, path="", **kwargs):
            return _client(request).storage.list_dir(
                f"{experiment_id}/{path}".rstrip("/")
            )

        @staticmethod
        def get_download_url(request, path=None, **kwargs):
            return f"/sdk/download?path={_path(path, **kwargs)}"

        @staticmethod
        def get_data_product_metadata(request, path=None, **kwargs):
            return {}

        @staticmethod
        def update_data_product_content(
            request,
            path=None,
            fileContentText=None,
            **kwargs,
        ):
            resolved_path = _path(path, **kwargs)
            content = fileContentText
            if content is None:
                content = kwargs.get("content")
            return _data_product(
                _client(request).storage.upload_file(
                    path=resolved_path,
                    content=content,
                    name=resolved_path.rsplit("/", 1)[-1],
                )
            )

    user_storage = _UserStorage()

    class _IntermediateOutput:
        @staticmethod
        def can_fetch_intermediate_output(request, experiment_model, output_name):
            return False

        @staticmethod
        def fetch_intermediate_output(request, experiment_id, output_name):
            return None

        @staticmethod
        def get_intermediate_output_data_products(
            request,
            experiment_model,
            output_name,
        ):
            return []

    class _ExperimentUtil:
        intermediate_output = _IntermediateOutput()

        @staticmethod
        def launch(request, experiment_id, **kwargs):
            from django.conf import settings

            return _client(request).research.launch_experiment(
                experiment_id,
                settings.GATEWAY_ID,
            )

    experiment_util = _ExperimentUtil()

    class _RemoteApi:
        @staticmethod
        def is_remote_api_configured():
            return False

        @staticmethod
        def call(*args, **kwargs):
            raise RuntimeError("remoteapi is not configured in the gRPC portal")

    remoteapi = _RemoteApi()
