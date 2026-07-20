"""Collect and verify scientific outputs for one portal run."""

import logging
import os
import re

from .airavata_grpc import user_storage
from . import (
    models,
    output_presentation_contracts,
    scientific_output_contracts,
    workflow_continuation as workflow_continuation_domain,
)


logger = logging.getLogger(__name__)

SCIENTIFIC_LOG_HEAD_BYTES = 256 * 1024
SCIENTIFIC_LOG_TAIL_BYTES = 768 * 1024
TERMINAL_SCHEDULER_STATUSES = {"COMPLETED", "FAILED", "CANCELED"}


def _is_epolyscat_run(run):
    applications = (
        run.module_application,
        run.workflow_application,
        run.utility_application,
    )
    return run.workflow_stage == "ePolyScat_Run" or any(
        re.sub(r"[^a-z0-9]", "", str(application or "").lower())
        == "epolyscat"
        for application in applications
    )


def epolyscat_output_declarations(request, run):
    if not _is_epolyscat_run(run):
        return {}

    declarations = {}
    script_inputs = models.Input.objects.filter(
        run=run,
        type="files",
        name__in=("ePolyscat_Input_File", "ePolyScat_Input_File"),
    ).prefetch_related("files")
    for input_instance in script_inputs:
        for file_instance in input_instance.files.all():
            file_object = None
            try:
                file_object = user_storage.open_file(
                    request,
                    data_product_uri=file_instance.data_product_uri,
                )
                declarations.update(
                    output_presentation_contracts.parse_file_name_declarations(
                        file_object.read()
                    )
                )
            except Exception:
                logger.debug(
                    "Unable to read ePolyScat output declarations from %s",
                    file_instance.name,
                    exc_info=True,
                )
            finally:
                close = getattr(file_object, "close", None)
                if callable(close):
                    close()
    return declarations


def output_files_for_run(request, run):
    output_files = []
    seen_file_names = set()
    seen_directories = set()

    def append_output_files(files):
        for file_data in files:
            if isinstance(file_data, dict):
                filename = (
                    file_data.get("name")
                    or file_data.get("filename")
                    or file_data.get("fileName")
                )
            else:
                filename = str(file_data)
            if filename and filename not in seen_file_names:
                output_files.append(file_data)
                seen_file_names.add(filename)

    def child_directory_path(parent_path, directory_data):
        if not isinstance(directory_data, dict):
            return ""
        directory_path = directory_data.get("path")
        if directory_path:
            return directory_path.strip("/")
        directory_name = (
            directory_data.get("name")
            or directory_data.get("filename")
            or directory_data.get("fileName")
        )
        if not directory_name and directory_data.get("resource_path"):
            directory_name = os.path.basename(
                directory_data["resource_path"].rstrip("/")
            )
        if not directory_name:
            return ""
        return "/".join(
            part.strip("/")
            for part in (parent_path, directory_name)
            if part and part.strip("/")
        )

    def append_experiment_directory_files(experiment_id, path="", depth=0):
        if depth > 4 or path in seen_directories:
            return
        seen_directories.add(path)
        try:
            directories, files = user_storage.list_experiment_dir(
                request,
                experiment_id,
                path=path,
            )
        except Exception:
            logger.debug(
                "Unable to list output directory %s for %s",
                path,
                experiment_id,
                exc_info=True,
            )
            return

        append_output_files(files)
        for directory_data in directories:
            child_path = child_directory_path(path, directory_data)
            if child_path:
                append_experiment_directory_files(
                    experiment_id,
                    path=child_path,
                    depth=depth + 1,
                )

    latest_execution = run.latest_execution
    if latest_execution is not None:
        experiment_id = latest_execution.airavata_experiment_id
        append_experiment_directory_files(experiment_id)
        append_experiment_directory_files(experiment_id, path="ARCHIVE")

    declared_file_types = epolyscat_output_declarations(request, run)
    return output_presentation_contracts.annotate_output_files(
        output_files,
        declared_file_types=declared_file_types,
    )


def _chunk_bytes(value):
    if isinstance(value, bytes):
        return value
    return str(value or "").encode("utf-8", errors="replace")


def _bounded_log_contents(file_object):
    head = _chunk_bytes(file_object.read(SCIENTIFIC_LOG_HEAD_BYTES))
    tail = bytearray()
    while True:
        chunk = file_object.read(64 * 1024)
        if not chunk:
            break
        tail.extend(_chunk_bytes(chunk))
        if len(tail) > SCIENTIFIC_LOG_TAIL_BYTES:
            del tail[: len(tail) - SCIENTIFIC_LOG_TAIL_BYTES]
    if not tail:
        contents = head
    else:
        contents = head + b"\n... [middle omitted] ...\n" + bytes(tail)
    return contents.decode("utf-8", errors="replace")


def read_manifest_entry_text(request, entry):
    data_product_uri = entry.get("data_product_uri", "")
    descriptor = entry.get("descriptor", {})
    path = ""
    if isinstance(descriptor, dict):
        path = descriptor.get("resource_path") or descriptor.get("path") or ""
    if not data_product_uri and not path:
        return None

    file_object = None
    try:
        if data_product_uri:
            file_object = user_storage.open_file(
                request,
                data_product_uri=data_product_uri,
            )
        else:
            file_object = user_storage.open_file(request, path=path)
        return _bounded_log_contents(file_object)
    except Exception:
        return None
    finally:
        close = getattr(file_object, "close", None)
        if callable(close):
            close()


def build_report(request, run, output_files=None):
    output_files = (
        output_files_for_run(request, run)
        if output_files is None
        else list(output_files)
    )
    input_files = run.inputs.filter(type="files").prefetch_related("files")
    input_file_names = {
        os.path.basename(file_instance.name)
        for input_instance in input_files
        for file_instance in input_instance.files.all()
    }
    input_data_product_uris = {
        file_instance.data_product_uri
        for input_instance in input_files
        for file_instance in input_instance.files.all()
        if file_instance.data_product_uri
    }

    def is_staged_input(file_data):
        if isinstance(file_data, dict):
            name = (
                file_data.get("name")
                or file_data.get("filename")
                or file_data.get("fileName")
                or ""
            )
            data_product_uri = (
                file_data.get("data_product_uri")
                or file_data.get("dataProductURI")
                or file_data.get("data-product-uri")
                or file_data.get("productUri")
                or file_data.get("product_uri")
                or file_data.get("value")
                or ""
            )
        else:
            name = str(file_data or "")
            data_product_uri = ""
        return (
            os.path.basename(name) in input_file_names
            or data_product_uri in input_data_product_uris
        )

    output_files = [
        file_data for file_data in output_files if not is_staged_input(file_data)
    ]
    manifest = scientific_output_contracts.build_output_manifest(output_files)
    classification = workflow_continuation_domain.classify_run(run) or {}
    source_application = classification.get("source_application") or (
        run.workflow_application
        or run.module_application
        or run.utility_application
    )
    content_cache = {}

    def read_text(entry):
        cache_key = entry.get("data_product_uri") or entry.get("name")
        if cache_key not in content_cache:
            content_cache[cache_key] = read_manifest_entry_text(request, entry)
        return content_cache[cache_key]

    verification = scientific_output_contracts.verify_scientific_outputs(
        source_application,
        manifest,
        read_text=read_text,
    )
    latest_execution = run.latest_execution
    scheduler_status = ""
    if latest_execution is not None:
        try:
            scheduler_status = str(
                latest_execution.get_airavata_experiment_status(request)
            ).upper()
        except Exception:
            scheduler_status = str(
                latest_execution.airavata_experiment_status or ""
            ).upper()
    return {
        "run_id": run.id,
        "source_stage": classification.get("source_stage", ""),
        "source_application": source_application,
        "scheduler_status": scheduler_status,
        "scheduler_complete": scheduler_status == "COMPLETED",
        "files": manifest,
        "scientific_verification": verification,
    }


def is_verified_complete(request, run):
    report = build_report(request, run)
    return bool(
        report["scheduler_complete"]
        and report["scientific_verification"]["status"] == "verified"
    )
