export function normalizeAnalysisApplications(applications, allowedApplications, fallbackApplication) {
  const allowed = new Set(allowedApplications || []);
  const normalized = [];

  (applications || []).forEach(applicationId => {
    if (allowed.has(applicationId) && !normalized.includes(applicationId)) {
      normalized.push(applicationId);
    }
  });

  if (normalized.length === 0 && allowed.has(fallbackApplication)) {
    normalized.push(fallbackApplication);
  }

  return normalized;
}

export function addAnalysisApplication(applications, applicationId, allowedApplications) {
  if (!(allowedApplications || []).includes(applicationId) || applications.includes(applicationId)) {
    return [...applications];
  }

  return [...applications, applicationId];
}

export function removeAnalysisApplication(applications, applicationId) {
  if (applications.length <= 1) {
    return [...applications];
  }

  return applications.filter(item => item !== applicationId);
}

export function moveAnalysisApplication(applications, index, offset) {
  const targetIndex = index + offset;
  if (index < 0 || index >= applications.length || targetIndex < 0 || targetIndex >= applications.length) {
    return [...applications];
  }

  const reordered = [...applications];
  const [applicationId] = reordered.splice(index, 1);
  reordered.splice(targetIndex, 0, applicationId);
  return reordered;
}
