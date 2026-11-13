// Local contracts describe supported inputs and stage ordering. Registration
// and application identity come exclusively from the discovered catalog.
export function filterApplicationContracts(contracts, discoveredApplications) {
  const discovered = new Map((discoveredApplications || []).map(app => [app.id, app]));
  return contracts.map(group => ({
    ...group,
    applications: group.id === "workflow"
      ? group.applications.map(stage => ({
        ...stage,
        workflowApplicationIds: (stage.workflowApplicationIds || []).filter(id => discovered.has(id)),
      }))
      : group.applications.filter(app => discovered.has(app.id) && discovered.get(app.id).kind === group.id)
        .map(app => ({ ...app, ...discovered.get(app.id) })),
  }));
}

export function applicationSelectionForRun(run) {
  const runType = run.runMode || "module";
  const metadata = run.workflowMetadata || {};
  if (runType === "workflow") {
    const stage = run.workflowStage || "Data_Gen";
    const applications = {
      Data_Gen: metadata.dataGenerationApplication || "",
      ePolyScat_Run: "ePolyScat",
      Analysis: (metadata.analysisApplications || [])[0] || "",
    };
    if (run.workflowStage) applications[stage] = run.workflowApplication || run.moduleApplication || run.utilityApplication || "";
    return { runType, applicationId: stage, stageSelections: applications, analysisApplications: metadata.analysisApplications || [] };
  }
  return { runType, applicationId: runType === "utility" ? run.utilityApplication || "" : run.moduleApplication || "ePolyScat" };
}
