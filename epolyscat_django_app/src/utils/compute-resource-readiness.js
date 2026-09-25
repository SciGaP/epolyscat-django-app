function deploymentValue(deployment, camelCaseName, snakeCaseName) {
  if (!deployment) return "";
  return deployment[camelCaseName] || deployment[snakeCaseName] || "";
}

function deploymentExecutableBasename(deployment) {
  const executablePath = String(
      deploymentValue(deployment, "executablePath", "executable_path")
  );
  return executablePath.split(/[\\/]/).pop().toLowerCase();
}

export function isUtilityDispatcherDeployment(deployment) {
  const executableBasename = deploymentExecutableBasename(deployment);
  return executableBasename.includes("controller")
    || (
      executableBasename.includes("epolyscat")
      && executableBasename.includes("wrapper")
    );
}

export function filterEligibleApplicationDeployments(
    deployments,
    executionKind = "module"
) {
  const deploymentsByHost = new Map();
  (deployments || []).forEach(deployment => {
    const computeHostId = deploymentValue(
        deployment,
        "computeHostId",
        "compute_host_id"
    );
    if (!computeHostId) return;
    if (!deploymentsByHost.has(computeHostId)) {
      deploymentsByHost.set(computeHostId, []);
    }
    deploymentsByHost.get(computeHostId).push(deployment);
  });

  return Array.from(deploymentsByHost.values()).reduce((eligible, hostDeployments) => {
    if (hostDeployments.length !== 1) return eligible;

    const deployment = hostDeployments[0];
    const appDeploymentId = deploymentValue(
        deployment,
        "appDeploymentId",
        "app_deployment_id"
    );
    const executablePath = deploymentValue(
        deployment,
        "executablePath",
        "executable_path"
    );
    if (!appDeploymentId || !executablePath) return eligible;
    if (executionKind === "utility" && !isUtilityDispatcherDeployment(deployment)) {
      return eligible;
    }

    eligible.push(deployment);
    return eligible;
  }, []);
}

export function buildComputeResourceReadiness({
  applicationModuleId,
  applicationLabel = "Application",
  groupResourceProfileId,
  computeResourceId,
  computeResourceLabel,
  eligibleDeployments = [],
  scopeLoaded = false,
  loading = false,
  error = null,
} = {}) {
  const label = applicationLabel || "Application";
  const result = (status, message, ready = false) => ({ status, message, ready });

  if (!applicationModuleId) {
    return result(
        "configuration-missing",
        `${label} deployment configuration is unavailable.`
    );
  }
  if (!groupResourceProfileId) {
    return result(
        "allocation-required",
        `Select an allocation to check ${label} deployments.`
    );
  }
  if (loading || !scopeLoaded) {
    return result("checking", `Checking ${label} deployments...`);
  }
  if (error) {
    return result(
        "load-error",
        `Unable to load ${label} deployments. Try again.`
    );
  }
  if (!eligibleDeployments.length) {
    return result(
        "no-deployment",
        `No deployment for ${label} is available in this allocation.`
    );
  }
  if (!computeResourceId) {
    return result(
        "selection-required",
        `Select a deployed compute resource for ${label}.`
    );
  }

  const selectedDeployment = eligibleDeployments.find(deployment => {
    return deploymentValue(deployment, "computeHostId", "compute_host_id")
      === computeResourceId;
  });
  if (!selectedDeployment) {
    return result(
        "invalid-selection",
        `The selected compute resource is not deployed for ${label}.`
    );
  }

  return result(
      "ready",
      `${label} is deployed on ${computeResourceLabel || computeResourceId}.`,
      true
  );
}
