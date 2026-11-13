import subprocess
import textwrap
from pathlib import Path


APP_DIR = Path(__file__).resolve().parents[1]


def _run_node_script(script):
    return subprocess.run(
        ["node", "-e", script],
        cwd=APP_DIR,
        text=True,
        capture_output=True,
        check=False,
    )


def test_deployment_filter_requires_one_runnable_deployment_per_host():
    script = textwrap.dedent(
        """
        const fs = require("fs");
        const path = require("path");
        const babel = require("@babel/core");

        const sourcePath = path.join(
          process.cwd(), "src", "utils", "compute-resource-readiness.js"
        );
        const source = fs.readFileSync(sourcePath, "utf8");
        const compiled = babel.transformSync(source, {
          plugins: ["@babel/plugin-transform-modules-commonjs"],
        }).code;
        const moduleObject = { exports: {} };
        new Function("module", "exports", compiled)(moduleObject, moduleObject.exports);

        const { filterEligibleApplicationDeployments } = moduleObject.exports;
        const deployments = [
          {
            appDeploymentId: "direct",
            computeHostId: "direct-host",
            executablePath: "/apps/ePolyScat/bin/ePolyScat",
          },
          {
            appDeploymentId: "controller",
            computeHostId: "controller-host",
            executablePath: "/apps/ePolyScat/epolyscat_controller.sh",
          },
          {
            appDeploymentId: "wrapper",
            computeHostId: "wrapper-host",
            executablePath: "/apps/ePolyScat/epolyscat_apps_wrapper.sh",
          },
          {
            appDeploymentId: "missing-path",
            computeHostId: "missing-path-host",
            executablePath: "",
          },
          {
            appDeploymentId: "duplicate-a",
            computeHostId: "duplicate-host",
            executablePath: "/apps/ePolyScat/bin/ePolyScat",
          },
          {
            appDeploymentId: "duplicate-b",
            computeHostId: "duplicate-host",
            executablePath: "/apps/ePolyScat/bin/ePolyScat-2",
          },
        ];

        const moduleHosts = filterEligibleApplicationDeployments(deployments, "module")
          .map(deployment => deployment.computeHostId);
        const utilityHosts = filterEligibleApplicationDeployments(deployments, "utility")
          .map(deployment => deployment.computeHostId);

        if (JSON.stringify(moduleHosts) !== JSON.stringify([
          "direct-host", "controller-host", "wrapper-host",
        ])) {
          throw new Error(`Unexpected module deployments: ${JSON.stringify(moduleHosts)}`);
        }
        if (JSON.stringify(utilityHosts) !== JSON.stringify([
          "controller-host", "wrapper-host",
        ])) {
          throw new Error(`Unexpected utility deployments: ${JSON.stringify(utilityHosts)}`);
        }
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr


def test_readiness_explains_configuration_loading_empty_invalid_and_ready_states():
    script = textwrap.dedent(
        """
        const fs = require("fs");
        const path = require("path");
        const babel = require("@babel/core");

        const sourcePath = path.join(
          process.cwd(), "src", "utils", "compute-resource-readiness.js"
        );
        const source = fs.readFileSync(sourcePath, "utf8");
        const compiled = babel.transformSync(source, {
          plugins: ["@babel/plugin-transform-modules-commonjs"],
        }).code;
        const moduleObject = { exports: {} };
        new Function("module", "exports", compiled)(moduleObject, moduleObject.exports);

        const { buildComputeResourceReadiness } = moduleObject.exports;
        const deployment = {
          appDeploymentId: "deployment",
          computeHostId: "frontera",
          executablePath: "/apps/ePolyScat/bin/ePolyScat",
        };
        const base = {
          applicationModuleId: "epolyscat-module",
          applicationLabel: "ePolyScat",
          groupResourceProfileId: "allocation",
          computeResourceId: "",
          computeResourceLabel: "",
          eligibleDeployments: [deployment],
          scopeLoaded: true,
          loading: false,
          error: null,
        };
        const cases = [
          [{ ...base, applicationModuleId: "" }, "configuration-missing"],
          [{ ...base, groupResourceProfileId: "" }, "allocation-required"],
          [{ ...base, loading: true, scopeLoaded: false }, "checking"],
          [{ ...base, error: new Error("network") }, "load-error"],
          [{ ...base, eligibleDeployments: [] }, "no-deployment"],
          [base, "selection-required"],
          [{ ...base, computeResourceId: "expanse" }, "invalid-selection"],
          [{
            ...base,
            computeResourceId: "frontera",
            computeResourceLabel: "Frontera",
          }, "ready"],
        ];

        cases.forEach(([options, expectedStatus]) => {
          const result = buildComputeResourceReadiness(options);
          if (result.status !== expectedStatus) {
            throw new Error(`${expectedStatus} became ${JSON.stringify(result)}`);
          }
          if (result.ready !== (expectedStatus === "ready")) {
            throw new Error(`Unexpected ready flag for ${expectedStatus}`);
          }
          if (!result.message.includes("ePolyScat")) {
            throw new Error(`Application label missing from ${JSON.stringify(result)}`);
          }
        });
        """
    )

    result = _run_node_script(script)

    assert result.returncode == 0, result.stderr
