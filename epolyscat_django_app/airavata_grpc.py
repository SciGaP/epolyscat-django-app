"""Compatibility layer for Thrift and gRPC Airavata portal runtimes.

The dev portal currently runs the gRPC Airavata SDK, while local and older
portal installs may still have the Thrift SDK. Keep the rest of this app using
the old small surface: enum objects with ``.name``/``.value``, ``user_storage``,
``experiment_util``, and ``remoteapi``.
"""

from types import SimpleNamespace


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
    from airavata.model.status.ttypes import ExperimentState
    from airavata.model.workspace.ttypes import Project
    from airavata_django_portal_sdk import experiment_util, remoteapi, user_storage

    def django_user(request):
        return request.user

else:
    ExperimentModel = _exp_pb2.ExperimentModel
    UserConfigurationDataModel = _exp_pb2.UserConfigurationDataModel
    ComputationalResourceSchedulingModel = (
        _sched_pb2.ComputationalResourceSchedulingModel
    )
    Project = _ws_pb2.Project

    class _EnumMember:
        def __init__(self, name, value):
            self.name = name
            self.value = value

        def __int__(self):
            return self.value

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
