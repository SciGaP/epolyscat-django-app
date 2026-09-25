"""Discover registered applications before applying the portal's workflow contracts."""

import json
import logging
import re

from epolyscat_django_app import runtime_audit

logger = logging.getLogger(__name__)

# Semantic capabilities, not IDs from a particular gateway.
SUPPORTED_MODULES = {
    "ePolyScat": "EPOLYSCAT_APPLICATION_ID",
    "Gaussian16": "GAUSSIAN16_APPLICATION_ID",
    "OpenMolcas": "OPENMOLCAS_APPLICATION_ID",
}


class DiscoveryUnavailable(Exception):
    pass


def load_registry(request, gateway_id):
    if hasattr(request, "_epolyscat_registry"):
        return request._epolyscat_registry
    try:
        client, token = request.airavata_client, request.authz_token
        registry = {
            "modules": client.getAccessibleAppModules(token, gateway_id) or [],
            "interfaces": client.getAllApplicationInterfaces(token, gateway_id) or [],
            "deployments": client.getAllApplicationDeployments(token, gateway_id) or [],
        }
    except Exception as error:
        logger.exception("Unable to discover Airavata applications")
        raise DiscoveryUnavailable("Unable to load applications from Airavata. Please try again.") from error
    request._epolyscat_registry = registry
    return registry


def _value(record, name, default=None):
    return record.get(name, default) if isinstance(record, dict) else getattr(record, name, default)


def _normalized_name(value):
    return re.sub(r"[^a-z0-9]", "", str(value or "").lower())


def _utility_options(interface):
    selector = next((item for item in (_value(interface, "applicationInputs", []) or [])
                     if _value(item, "name") == "Application_Utility"), None)
    metadata = _value(selector, "metaData", {}) or {}
    if isinstance(metadata, str):
        try:
            metadata = json.loads(metadata)
        except (ValueError, TypeError):
            return set()
    if not isinstance(metadata, dict) or not isinstance(metadata.get("editor"), dict):
        return set()
    # Airavata's input editors use editor.config.options; older registrations
    # also stored options directly under editor.
    config = metadata["editor"].get("config", metadata["editor"])
    if not isinstance(config, dict):
        return set()
    options = config.get("options", [])
    if not isinstance(options, list):
        return set()
    return {str(option.get("value", "") if isinstance(option, dict) else option) for option in options}


def build_catalog(*, modules, interfaces, deployments, configured_ids=None):
    configured_ids = configured_ids or {}
    catalog = {"EPOLYSCAT": {key: None for key in SUPPORTED_MODULES.values()}, "applications": [], "issues": []}
    for capability, setting_key in SUPPORTED_MODULES.items():
        configured_id = configured_ids.get(setting_key)
        matches = [module for module in modules if (
            _value(module, "appModuleId") == configured_id if configured_id
            else _normalized_name(_value(module, "appModuleName")) == _normalized_name(capability)
        )]
        if len(matches) != 1:
            if configured_id or matches:
                catalog["issues"].append(f"{capability}: the configured application is unavailable or its registration is ambiguous.")
            continue
        module = matches[0]
        module_id = _value(module, "appModuleId")
        app_interfaces = [interface for interface in interfaces
                          if module_id in (_value(interface, "applicationModules", []) or [])]
        app_deployments = [deployment for deployment in deployments
                           if _value(deployment, "appModuleId") == module_id
                           and _value(deployment, "appDeploymentId")
                           and _value(deployment, "computeHostId")
                           and _value(deployment, "executablePath")]
        if not module_id or len(app_interfaces) != 1 or not _value(app_interfaces[0], "applicationInterfaceId") or not app_deployments:
            catalog["issues"].append(f"{capability}: one application interface and a deployment are required.")
            continue
        interface = app_interfaces[0]
        application = dict(id=capability, label=_value(module, "appModuleName") or capability,
                           kind="module", moduleId=module_id,
                           interfaceId=_value(interface, "applicationInterfaceId"))
        catalog["EPOLYSCAT"][setting_key] = module_id
        catalog["applications"].append(application)
        if capability == "ePolyScat":
            advertised = _utility_options(interface)
            audit = runtime_audit.audit_runtime_configuration(
                application_module_id=module_id, modules=modules,
                interfaces=interfaces, deployments=app_deployments,
            )
            catalog["applications"].extend(
                dict(application, id=utility["id"], label=utility["id"], kind="utility")
                for utility in audit["utilities"] if utility["ready"] and utility["id"] in advertised
            )
            for utility in audit["utilities"]:
                if utility["id"] not in advertised or utility["ready"]:
                    continue
                reasons = []
                if not utility["deployed"]:
                    reasons.append("no utility dispatcher is deployed")
                if not utility["control_input_configured"]:
                    reasons.append("a control-file input is missing")
                if not utility["required_inputs_configured"]:
                    reasons.append("required data inputs are missing")
                elif not utility["supports_required_file_count"]:
                    reasons.append("a multi-file input is not configured")
                catalog["issues"].append(f"{utility['id']}: " + "; ".join(reasons) + ".")
    return catalog
