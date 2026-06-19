"""gRPC compatibility layer for the new (gRPC) Airavata portal.

Re-exposes, against ``request.airavata`` (the gRPC AiravataClient), the surface
epolyscat used from the old Thrift SDK + airavata-django-portal-sdk: proto type
aliases, an ``ExperimentState`` with the old bare names, and ``user_storage`` /
``experiment_util`` shims.
"""
from airavata_sdk.generated.org.apache.airavata.model.experiment import experiment_pb2 as _exp_pb2
from airavata_sdk.generated.org.apache.airavata.model.scheduling import scheduling_pb2 as _sched_pb2
from airavata_sdk.generated.org.apache.airavata.model.status import status_pb2 as _status_pb2
from airavata_sdk.generated.org.apache.airavata.model.application.io import application_io_pb2 as _io_pb2
from airavata_sdk.generated.org.apache.airavata.model.workspace import workspace_pb2 as _ws_pb2

ExperimentModel = _exp_pb2.ExperimentModel
UserConfigurationDataModel = _exp_pb2.UserConfigurationDataModel
ComputationalResourceSchedulingModel = _sched_pb2.ComputationalResourceSchedulingModel
DataType = _io_pb2.DataType
Project = _ws_pb2.Project


class ExperimentState:
    """Old bare ``ExperimentState.<NAME>`` constants -> new proto values."""
    _E = _status_pb2.ExperimentState
    CREATED = _E.EXPERIMENT_STATE_CREATED
    VALIDATED = _E.EXPERIMENT_STATE_VALIDATED
    SCHEDULED = _E.EXPERIMENT_STATE_SCHEDULED
    LAUNCHED = _E.EXPERIMENT_STATE_LAUNCHED
    EXECUTING = _E.EXPERIMENT_STATE_EXECUTING
    CANCELING = _E.EXPERIMENT_STATE_CANCELING
    CANCELED = _E.EXPERIMENT_STATE_CANCELED
    COMPLETED = _E.EXPERIMENT_STATE_COMPLETED
    FAILED = _E.EXPERIMENT_STATE_FAILED

    _NAMES_TO_VALUES = {
        "CREATED": CREATED, "VALIDATED": VALIDATED, "SCHEDULED": SCHEDULED,
        "LAUNCHED": LAUNCHED, "EXECUTING": EXECUTING, "CANCELING": CANCELING,
        "CANCELED": CANCELED, "COMPLETED": COMPLETED, "FAILED": FAILED,
    }
    _VALUES_TO_NAMES = {v: k for k, v in _NAMES_TO_VALUES.items()}


def _client(request):
    return request.airavata


def django_user(request):
    """Bridge the Keycloak identity (``request.user``, no DB row) to a Django
    ``auth.User``, which epolyscat's models use as the ``owner`` FK."""
    from django.contrib.auth import get_user_model
    u = getattr(request, "user", None)
    username = getattr(u, "username", None)
    if not username:
        return u
    User = get_user_model()
    obj, _ = User.objects.get_or_create(
        username=username,
        defaults={
            "email": getattr(u, "email", "") or "",
            "first_name": getattr(u, "first_name", "") or "",
            "last_name": getattr(u, "last_name", "") or "",
        },
    )
    return obj


class _UserStorage:
    """Path-based, backed by ``request.airavata.storage``."""

    @staticmethod
    def exists(request, path, **kw):
        return _client(request).storage.file_exists(path)

    @staticmethod
    def user_file_exists(request, path, **kw):
        try:
            return _client(request).storage.file_exists(path)
        except Exception:
            return False

    @staticmethod
    def dir_exists(request, path, **kw):
        return _client(request).storage.dir_exists(path)

    @staticmethod
    def open_file(request, path, **kw):
        import io
        return io.BytesIO(_client(request).storage.download_file(path))

    @staticmethod
    def save(request, path, file, name=None, **kw):
        content = file.read() if hasattr(file, "read") else file
        name = name or (getattr(file, "name", None) or path.rsplit("/", 1)[-1])
        return _client(request).storage.upload_file(path=path, content=content, name=name)

    @staticmethod
    def save_input_file(request, file, name=None, **kw):
        content = file.read() if hasattr(file, "read") else file
        name = name or getattr(file, "name", "input")
        return _client(request).storage.upload_file(path=name, content=content, name=name)

    @staticmethod
    def create_user_dir(request, path, **kw):
        return _client(request).storage.create_dir(path)

    @staticmethod
    def delete(request, path, **kw):
        return _client(request).storage.delete_file(path)

    @staticmethod
    def delete_dir(request, path, **kw):
        return _client(request).storage.delete_dir(path)

    @staticmethod
    def listdir(request, path, **kw):
        return _client(request).storage.list_dir(path)

    @staticmethod
    def list_experiment_dir(request, experiment_id, path="", **kw):
        return _client(request).storage.list_dir(f"{experiment_id}/{path}".rstrip("/"))

    @staticmethod
    def get_download_url(request, path, **kw):
        return f"/sdk/download?path={path}"

    @staticmethod
    def get_data_product_metadata(request, path, **kw):
        return {}

    @staticmethod
    def update_data_product_content(request, path, content, **kw):
        name = path.rsplit("/", 1)[-1]
        return _client(request).storage.upload_file(path=path, content=content, name=name)


user_storage = _UserStorage()


class _ExperimentUtil:
    """Backed by ``request.airavata.research``."""

    @staticmethod
    def launch(request, experiment_id, **kw):
        from django.conf import settings
        return _client(request).research.launch_experiment(experiment_id, settings.GATEWAY_ID)

    @staticmethod
    def intermediate_output(request, experiment_id, output_names=None, **kw):
        return _client(request).research.get_intermediate_outputs(experiment_id, output_names)


experiment_util = _ExperimentUtil()
