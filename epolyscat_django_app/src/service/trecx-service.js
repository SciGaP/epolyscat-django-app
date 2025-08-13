import axios from "axios";
import http from "http";
import https from "https";

const httpAgent = new http.Agent({keepAlive: true});
const httpsAgent = new https.Agent({keepAlive: true});

const {utils} = AiravataAPI;

const axiosInstance = axios.create({
    httpAgent,
    httpsAgent,
    baseURL: "/",
    withCredentials: false,
    headers: {
        'Accept': '*/*',
        'Content-Type': 'application/json',
        'Access-Control-Allow-Origin': "*",
        'Origin': '*',
        'x-csrftoken': utils.FetchUtils.getCSRFToken()
    }
})

export const ExperimentService = {
    encodeObj(obj) {
        return {
            experimentId: obj.id,
            name: obj.name,
            description: obj.description,
            owner: obj.owner,
            created: new Date(obj.created).toLocaleString(),
            updated: new Date(obj.updated).toLocaleString(),
            deleted: obj.deleted,
            airavataProjectId: obj.airavata_project_id,
            activeRunCount: obj.active_run_count,
            runCount: obj.run_count,
            root: obj.root,
        };
    },
    async fetchExperimentStatistics() {
        const {data: {experiments_count, runs_count}} = await axiosInstance.get(`/trecx_django_app/api/experiments/statistics/`);
        return {
            experimentCount: experiments_count,
            runCount: runs_count
        }
    },
    async fetchAllExperiments({page = 1, pageSize = 1000} = {page: 1, pageSize: 1000}) {
        const {data} = await axiosInstance.get(`/trecx_django_app/api/experiments/?page=${page}&page_size=${pageSize}`);
        data.results = data.results.map(this.encodeObj);
        return data;
    },
    async fetchExperiment({experimentId}) {
        const {data} = await axiosInstance.get(`/trecx_django_app/api/experiments/${experimentId}`);
        return this.encodeObj(data);
    },
    async createExperiment({name, description}) {
        let {data} = await axiosInstance.post("/trecx_django_app/api/experiments/", {
            "name": name,
            "description": description
        });

        return this.encodeObj(data);
    },
    async deleteExperiment({experimentId = null} = {}) {
        await axiosInstance.delete(`/trecx_django_app/api/experiments/${experimentId}/`);
    }
}

export const PlotService = {
    async plotSelectedRuns({runIds, expectationValue, xAxis, yAxis, flags}) {
        const {data} = await axiosInstance.post("/trecx_django_app/api/plot/", {
            "runs": runIds,
            "plotfile": expectationValue,
            "plot_parameters": {"xaxis": `${xAxis}`, "yaxes": `${yAxis}`, "flags": `${flags}`}
        });

        return {
            "plotImageUrl": data["mime-type"] && data["plot"] ? `data:${data["mime-type"]};base64,${data["plot"]}` : null,
            "output": data["output"],
            "userGuidance": data["user_guidance"]
        }
    },
    async getRunListInputs({runIds}) {
        const {data} = await axiosInstance.post("/trecx_django_app/api/list-inputs/", {
            "runs": runIds
        });

        return {
            "output": data.output
        };
    },
    async getRunListInputDifference({runIds}) {
        const {data} = await axiosInstance.post("/trecx_django_app/api/diff-inputs/", {
            "runs": runIds
        });

        return {
            "output": data.output
        };
    },
    async getPlotables({runIds}) {
        const {data} = await axiosInstance.post("/trecx_django_app/api/plotables/", {
            "runs": runIds
        })

        return data.filenames;
    },
    async getPlotParameters() {
        let {data} = await axiosInstance.get("/trecx_django_app/api/plot-parameters/");
        data = data.map(({xaxis, yaxes, flags}) => {
            return {text: `x=${xaxis} y=${yaxes} ${flags}`, value: {xAxis: xaxis, yAxis: yaxes, flags: flags}};
        });

        return data;
    },
    async getViewables({runId}) {
        const {data} = await axiosInstance.get(`/trecx_django_app/api/runs/${runId}/viewables/`)

        return data;
    },
    async getInputFiles({runId}) {
        const {data} = await axiosInstance.get(`/trecx_django_app/api/runs/${runId}/input-files/`)

        return data;
    },
}

export const RunService = {
    submitAllowedStatuses: ["FAILED", "Unsubmitted"],

    encodeObj(obj) {
        return {
            runId: obj.id,
            name: obj.name,
            number: obj.number,
            root: obj.root,
            experimentId: obj.experiment,
            created: new Date(obj.created).toLocaleString(),
            updated: new Date(obj.updated).toLocaleString(),
            deleted: obj.deleted,
            inpcDownloadUrl: obj.inpc_download_url,
            groupResourceProfileId: obj.group_resource_profile_id,
            computeResourceId: obj.compute_resource_id,
            queueName: obj.queue_name,
            coreCount: obj.core_count,
            nodeCount: obj.node_count,
            wallTimeLimit: obj.walltime_limit,
            totalPhysicalMemory: obj.total_physical_memory,
            status: obj.status,
            resource: obj.resource,
            resourceShort: obj.resource_short,
            executions: obj.executions,
            inputTable: obj.input_table,
            canResubmit: obj.can_resubmit,
            canSubmit: RunService.submitAllowedStatuses.indexOf(obj.status) >= 0
        };
    },
    async fetchViewableContent({runId, filename, inpcDownloadUrl}) {
        let res;
        if (filename === "inpc" && inpcDownloadUrl) {
            res = await axios.get(inpcDownloadUrl);
        } else {
            res = await axiosInstance.get(`/trecx_django_app/api/runs/${runId}/viewables/${filename}/`);
        }

        return res.data;
    },
    async fetchAllRuns({experimentId = null, viewId = null, page = 1, pageSize = 1000} = {page: 1, pageSize: 1000}) {
        let queryString = `?page=${page}&page_size=${pageSize}`;
        if (experimentId) {
            queryString += `&experiment=${experimentId}`
        }

        if (viewId) {
            queryString += `&viewId=${viewId}`
        }

        const {data} = await axiosInstance.get("/trecx_django_app/api/runs/" + queryString);
        data.results = data.results.map(this.encodeObj);

        return data;
    },
    async fetchRun({runId = null} = {}) {
        const {data} = await axiosInstance.get(`/trecx_django_app/api/runs/${runId}`);
        const _run = this.encodeObj(data);

        if (!_run.inputTable) {
            const _clone = await this.cloneRun({runId});
            _run.inputTable = _clone.inputTable;
        }

        return _run;
    },
    async cloneRun({runId = null} = {}) {
        const {data} = await axiosInstance.post(`/trecx_django_app/api/runs/${runId}/new/`);

        return this.encodeObj(data);
    },
    async deleteRun({runId = null} = {}) {
        await axiosInstance.delete(`/trecx_django_app/api/runs/${runId}/`);
    },
    async createRun({root, experimentId, directedit, inputTable, groupResourceProfileId, computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory}, submit = false) {
        let data = {
            "root": root,
            "experiment": experimentId,
            "directedit": directedit,
            "input_table": inputTable,
            "group_resource_profile_id": groupResourceProfileId,
            "compute_resource_id": computeResourceId,
            "queue_name": queueName,
            "core_count": coreCount,
            "node_count": nodeCount,
            "walltime_limit": wallTimeLimit,
            "total_physical_memory": totalPhysicalMemory
        }

        const runCreateRes = await axiosInstance.post("/trecx_django_app/api/runs/", data);
        data = {...data, ...runCreateRes.data};

        if (submit) {
            const runSubmitRes = await axiosInstance.post(`/trecx_django_app/api/runs/${data.id}/submit/`, data);
            data = {...data, ...runSubmitRes.data};
        }

        return this.encodeObj(data);
    },
    async submitRun({runId, groupResourceProfileId, computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory}) {
        let data = {
            "group_resource_profile_id": groupResourceProfileId,
            "compute_resource_id": computeResourceId,
            "queue_name": queueName,
            "core_count": coreCount,
            "node_count": nodeCount,
            "walltime_limit": wallTimeLimit,
            "total_physical_memory": totalPhysicalMemory
        }

        const runSubmitRes = await axiosInstance.post(`/trecx_django_app/api/runs/${runId}/submit/`, data);
        data = {...data, ...runSubmitRes.data};

        return this.encodeObj(data);
    },
    async resubmitRun({runId, groupResourceProfileId, computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory}) {
        let data = {
            "group_resource_profile_id": groupResourceProfileId,
            "compute_resource_id": computeResourceId,
            "queue_name": queueName,
            "core_count": coreCount,
            "node_count": nodeCount,
            "walltime_limit": wallTimeLimit,
            "total_physical_memory": totalPhysicalMemory
        }

        const runResubmitRes = await axiosInstance.post(`/trecx_django_app/api/runs/${runId}/resubmit/`, data);

        data = {...data, ...runResubmitRes.data};

        return this.encodeObj(data);
    },
}

export const SettingsService = {
    async all() {
        const {data} = await axiosInstance.get("/trecx_django_app/api/settings/");
        return data;
    },
};

export const ViewService = {
    readonlyViewTypes: ["unsubmitted", "tutorial"],

    encodeObj(obj) {
        return {
            viewId: obj.id,
            name: obj.name,
            owner: obj.owner,
            created: new Date(obj.created).toLocaleString(),
            updated: new Date(obj.updated).toLocaleString(),
            deleted: obj.deleted,
            type: obj.type,
            activeRunCount: obj.active_run_count,
            runCount: obj.run_count,
            readonly: ViewService.readonlyViewTypes.indexOf(obj.type) >= 0
        };
    },
    async fetchAllViews({page = 1, pageSize = 1000, tutorials = false} = {page: 1, pageSize: 1000, tutorials: false}) {
        let url, data, res;
        if (tutorials) {
            url = "/epolyscat_django_app/api/views/tutorials/";
            res = await axiosInstance.get(url);
            data = {"count": 1, "next": null, "previous": null, "results": [res.data].map(this.encodeObj)}
        } else {
            url = `/trecx_django_app/api/views/?page=${page}&page_size=${pageSize}`;
            res = await axiosInstance.get(url);
            data = res.data;
            data.results = data.results.map(this.encodeObj);
        }

        return data;
    },
    async fetchView({viewId}) {
        const {data} = await axiosInstance.get(`/trecx_django_app/api/views/${viewId}`);
        return this.encodeObj(data);
    },
    async createView({name, runIds}) {
        let {data} = await axiosInstance.post("/trecx_django_app/api/views/", {
            "name": name
        });
        const {viewId} = this.encodeObj(data);
        await axiosInstance.put(`/trecx_django_app/api/views/${viewId}/add-runs/`, {
            "runs": runIds
        });

        return this.encodeObj({id: viewId, name, runs: runIds});
    },
    async updateView({viewId, name, runIds}) {
        await axiosInstance.put(`/trecx_django_app/api/views/${viewId}/`, {
            "name": name
        });
        await axiosInstance.put(`/trecx_django_app/api/views/${viewId}/add-runs/`, {
            "runs": runIds
        });

        return this.encodeObj({id: viewId, name, runs: runIds});
    },
    async removeRuns({viewId, runIds}) {
        await axiosInstance.put(`/trecx_django_app/api/views/${viewId}/remove-runs/`, {
            "runs": runIds
        });
    },
    async addRuns({viewId, runIds}) {
        await axiosInstance.put(`/trecx_django_app/api/views/${viewId}/add-runs/`, {
            "runs": runIds
        });
    },
    async deleteView({viewId = null} = {}) {
        await axiosInstance.delete(`/trecx_django_app/api/views/${viewId}/`);
    }
}
