import { fileMetadata } from "@/fileData";
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

const appBaseUrl = "/epolyscat_django_app/api/";

export const InputService = {
    async fetchApplicationInputs() {
        // const { data: { applicationInputs } } = await axiosInstance.get("/api/applications/ePolyScat_940ab1c9-4ceb-431c-8595-c6246a195442/application_interface/")
        
        return [
            {
                "name": "Calculation_Type",
                "value": null,
                "type": 0,
                "applicationArgument": null,
                "standardInput": false,
                "userFriendlyDescription": null,
                "metaData": {
                    "editor": {
                        "ui-component-id": "radio-button-input-editor",
                        "config": {
                            "options": [
                                {
                                    "value": "MODULE",
                                    "text": "MODULE"
                                },
                                {
                                    "value": "UTILITY",
                                    "text": "UTILITY"
                                },
                                {
                                    "value": "WORKFLOW",
                                    "text": "WORKFLOW"
                                }
                            ]
                        }
                    }
                },
                "inputOrder": 0,
                "isRequired": true,
                "requiredToAddedToCommandLine": true,
                "dataStaged": false,
                "storageResourceId": null,
                "isReadOnly": false,
                "overrideFilename": null
            },
            {
                "name": "Application_Workflow",
                "value": "",
                "type": 0,
                "applicationArgument": null,
                "standardInput": false,
                "userFriendlyDescription": "",
                "metaData": {
                    "editor": {
                        "ui-component-id": "radio-button-input-editor",
                        "config": {
                            "options": [
                                {
                                    "value": "Data_Gen",
                                    "text": "Data_Generation"
                                },
                                {
                                    "value": "ePolyScat_Run",
                                    "text": "ePolyScat_Run"
                                },
                                {
                                    "value": "Analysis",
                                    "text": "Analysis"
                                }
                            ]
                        },
                        "dependencies": {
                            "show": {
                                "OR": [
                                    {
                                        "Calculation_Type": {
                                            "comparison": "equals",
                                            "value": "WORKFLOW"
                                        }
                                    }
                                ]
                            },
                            "showOptions": {
                                "isRequired": true
                            }
                        }
                    }
                },
                "inputOrder": 1,
                "isRequired": true,
                "requiredToAddedToCommandLine": true,
                "dataStaged": false,
                "storageResourceId": null,
                "isReadOnly": false,
                "overrideFilename": null
            },
            {
                "name": "EPOLYSCAT_Application_Module",
                "value": null,
                "type": 0,
                "applicationArgument": null,
                "standardInput": false,
                "userFriendlyDescription": null,
                "metaData": {
                    "editor": {
                        "ui-component-id": "radio-button-input-editor",
                        "config": {
                            "options": [
                                {
                                    "value": "Gaussian16",
                                    "text": "Gaussian16"
                                },
                                {
                                    "value": "ePolyScat",
                                    "text": "ePolyScat"
                                },
                                {
                                    "value": "OpenMolcas",
                                    "text": "OpenMolcas"
                                }
                            ]
                        },
                        "dependencies": {
                            "show": {
                                "OR": [
                                    {
                                        "Calculation_Type": {
                                            "comparison": "equals",
                                            "value": "MODULE"
                                        }
                                    }
                                ]
                            },
                            "showOptions": {
                                "isRequired": true
                            }
                        }
                    }
                },
                "inputOrder": 2,
                "isRequired": true,
                "requiredToAddedToCommandLine": true,
                "dataStaged": false,
                "storageResourceId": null,
                "isReadOnly": false,
                "overrideFilename": null
            },
/*
            {
                "name": "Data_Gen",
                "value": null,
                "type": 0,
                "applicationArgument": null,
                "standardInput": false,
                "userFriendlyDescription": null,
                "metaData": {
                    "editor": {
                        "ui-component-id": "radio-button-input-editor",
                        "config": {
                            "options": [
                                {
                                    "value": "Gaussian16",
                                    "text": "Gaussian16"
                                },
                                {
                                    "value": "OpenMolcas",
                                    "text": "OpenMolcas"
                                }
                            ]
                        },
                        "dependencies": {
                                "show": {
                                    "Application_Workflow": {
                                        "comparison": "equals",
                                        "value": "Data_Gen"
                                    }
                                },
                                "showOptions": {
                                    "isRequired": true
                              }
                        },

                    }
                },
                "inputOrder": 3,
                "isRequired": true,
                "requiredToAddedToCommandLine": true,
                "dataStaged": false,
                "storageResourceId": null,
                "isReadOnly": false,
                "overrideFilename": null
            },
*/
            {
                "name": "Application_Utility",
                "value": null,
                "type": 0,
                "applicationArgument": null,
                "standardInput": false,
                "userFriendlyDescription": null,
                "metaData": {
                    "editor": {
                        "ui-component-id": "radio-button-input-editor",
                        "config": {
                            "options": [
                                {
                                    "value": "CnvMath",
                                    "text": "ConvertToMathematica"
                                },
                                {
                                    "value": "CnvMatLab",
                                    "text": "ConvertToMatlab"
                                },
                                {
                                    "value": "CnvLinFull",
                                    "text": "Compute Ph.ion.diff.Xsec."
                                },
                                {
                                    "value": "MoldenMerge",
                                    "text": "Merge Molde Data files"
                                },
                                {
                                    "value": "NRFPAD",
                                    "text": "computes the N photon RFPAD from the diff. Xsection"
                                },
                                {
                                    "value": "Cube2igor",
                                    "text": "G16CubeToIGOR_Plots"
                                }
                            ]
                        },
                        "dependencies": {
                            "show": {
                                "OR": [
                                    {
                                        "Calculation_Type": {
                                            "comparison": "equals",
                                            "value": "UTILITY"
                                        }
                                    }
                                ]
                            },
                            "showOptions": {
                                "isRequired": true
                            }
                        }
                    }
                },
                "inputOrder": 4,
                "isRequired": true,
                "requiredToAddedToCommandLine": true,
                "dataStaged": false,
                "storageResourceId": null,
                "isReadOnly": false,
                "overrideFilename": null
            },
            {
                "name": "Gaussian_Input",
                "value": null,
                "type": 3,
                "applicationArgument": null,
                "standardInput": true,
                "userFriendlyDescription": "File containing Gaussian calculation parameters.",
                "metaData": {
                    "editor": {
                        "dependencies": {
                            "show": {
                                "OR": [
                                    {
                                        "EPOLYSCAT_Application_Module": {
                                            "comparison": "equals",
                                            "value": "Gaussian16"
                                        }
                                        
                                    }
                                    /* ,                                                      
                                    {
                                        "Data_Gen": {
                                            "comparison": "equals",
                                            "value": "Gaussian16"
                                        }                                       
                                    } */                                   
                                ]                                                                        
                            },
                            "showOptions": {
                                "isRequired": true
                            }
                        }                    
                    }
                },
                "inputOrder": 5,
                "isRequired": true,
                "requiredToAddedToCommandLine": false,
                "dataStaged": true,
                "storageResourceId": null,
                "isReadOnly": false,
                "overrideFilename": "gaussian.inp"
            },
            {
                "name": "Molcas_Input",
                "value": "",
                "type": 3,
                "applicationArgument": null,
                "standardInput": true,
                "userFriendlyDescription": "Molcas Input file",
                "metaData": {
                    "editor": {
                        "dependencies": {
                            "show": {
                                "OR": [
                                    {
                                        "EPOLYSCAT_Application_Module": {
                                            "comparison": "equals",
                                            "value": "OpenMolcas"
                                        }
                                        
                                    }
                                    /*,                                    
                                    {
                                        "Data_Gen": {
                                            "comparison": "equals",
                                            "value": "OpenMolcas"
                                        }
                                    } */                                                                                                           
                                ]
                            },
                            "showOptions": {
                                "isRequired": true
                            }
                        }
                    }
                },
                "inputOrder": 6,
                "isRequired": true,
                "requiredToAddedToCommandLine": false,
                "dataStaged": true,
                "storageResourceId": null,
                "isReadOnly": false,
                "overrideFilename": null
            },
            {
                "name": "ePolyScat_Input_Data",
                "value": null,
                "type": 3,
                "applicationArgument": null,
                "standardInput": true,
                "userFriendlyDescription": "ePolyScat Input Data from programs such as OpenMolcas/Gaussian16/Molden.",
                "metaData": {
                    "editor": {
                        "dependencies": {
                            "show": {
                                "OR": [
                                    {
                                        "EPOLYSCAT_Application_Module": {
                                            "comparison": "equals",
                                            "value": "ePolyScat"
                                        }
                                    },
                                    {
                                        "Application_Workflow": {
                                            "comparison": "equals",
                                            "value": "ePolyScat_Run"
                                        }
                                    }
                                ]
                            },
                            "showOptions": {
                                "isRequired": true
                            }
                        }
                    }
                },
                "inputOrder": 7,
                "isRequired": true,
                "requiredToAddedToCommandLine": true,
                "dataStaged": true,
                "storageResourceId": null,
                "isReadOnly": false,
                "overrideFilename": null
            },
            {
                "name": "ePolyscat_Input_File",
                "value": "",
                "type": 3,
                "applicationArgument": null,
                "standardInput": true,
                "userFriendlyDescription": "This file contains the commands for the ePolyScat program",
                "metaData": {
                    "editor": {
                        "dependencies": {
                            "show": {
                                "OR": [
                                    {
                                        "EPOLYSCAT_Application_Module": {
                                            "comparison": "equals",
                                            "value": "ePolyScat"
                                        }
                                    },
                                    {
                                        "Application_Workflow": {
                                            "comparison": "equals",
                                            "value": "ePolyScat_Run"
                                        }
                                    }
                                ]
                            },
                            "showOptions": {
                                "isRequired": true
                            }
                        }
                    }
                },
                "inputOrder": 9,
                "isRequired": true,
                "requiredToAddedToCommandLine": true,
                "dataStaged": true,
                "storageResourceId": null,
                "isReadOnly": false,
                "overrideFilename": null
            },
            {
                "name": "molden.dat",
                "value": "",
                "type": 3,
                "applicationArgument": null,
                "standardInput": true,
                "userFriendlyDescription": "Molden Data file containing wavefunction.",
                "metaData": {
                    "editor": {
                        "dependencies": {
                            "show": {
                                "OR": [
                                    /*
                                    {
                                        "EPOLYSCAT_Application_Module": {
                                            "comparison": "equals",
                                            "value": "ePolyScat"
                                        }
                                    },
                                    {
                                        "Application_Workflow": {
                                            "comparison": "equals",
                                            "value": "ePolyScat_Run"
                                        }
                                    },
                                    */
                                    {
                                        "Application_Utility": {
                                            "comparison": "equals",
                                            "value": "MoldenMerge"
                                        }
                                    }
                                ]
                            },
                            "showOptions": {
                                "isRequired": true
                            }
                        }
                    }
                },
                "inputOrder": 8,
                "isRequired": true,
                "requiredToAddedToCommandLine": false,
                "dataStaged": true,
                "storageResourceId": null,
                "isReadOnly": false,
                "overrideFilename": null
            }
            
        ]
        
        // return applicationInputs;
    },
    async fetchOutputs(runId) {
        const { data } = await axiosInstance.get(appBaseUrl + `runs/${runId}/get_output_files/`, {});

        for (const outputFileObj of data)
            outputFileObj.downloadURL = `sdk/download/?data-product-uri=${encodeURIComponent(outputFileObj["data-product-uri"])}`;

        return data;
    },
    async fetchFileContents(file) {
        if (file == null)
            return null;

        if ("contents" in file)
            return file.contents;

        if ("base64Contents" in file)
            return null;

        if (file.isFromComputer)
            return await new Promise((resolve, reject) => {
                let reader = new FileReader();

                reader.onload = () => resolve(this.getValidatedResult(file.name, reader.result));

                reader.onerror = reject;

                reader.readAsText(file);
            });

        if (!("mimeType" in file) || (file["mimeType"] || "").startsWith("text")) {
            return this.getValidatedResult(file.name, (await axiosInstance.get(file.downloadURL)).data);
        }

        return null;
    },
    getValidatedResult(filename, fileContents) {
        const isPlaintext = fileMetadata.isPlaintext(filename);

        if (isPlaintext === false || (isPlaintext == null && /\0|�/.test(fileContents)))
            return null;
        else
            return fileContents;
    },
    async fetchFileContentsBase64(file) {
        if ("contents" in file)
            return btoa(file.contents);

        let reader = new FileReader();

        if (file.isFromComputer) {
            reader.readAsDataURL(file);
        } else {
            const { data } = await axiosInstance.get(file.downloadURL, { responseType: 'arraybuffer' });

            reader.readAsDataURL(new Blob([data]));
        }

        return await new Promise((resolve, reject) => {
            reader.onload = () => {
                resolve(reader.result.split("base64,")[1])
            };

            reader.onerror = reject;
        });
    },
    getProperFilenameOf(filename, files) {
        if (files.indexOf(filename) >= 0)
            return filename
        else if (filename == "target.0")
            return "target"

        const nnnFilename = filename.replaceAll(/\d{3}/g, "nnn");
        const isCFile = filename.slice(-2) == ".c";
        const isBswFile = filename.slice(-4) == ".bsw";

        if (files.indexOf(nnnFilename) >= 0)
            return nnnFilename;
        else if (isCFile && files.indexOf("name.c") >= 0)
            return "name.c"
        else if (isBswFile && files.indexOf("name.bsw") >= 0)
            return "name.bsw"

        return null;
    }
}

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
        const {data: {experiments_count, runs_count}} = await axiosInstance.get(`/epolyscat_django_app/api/experiments/statistics/`);
        return {
            experimentCount: experiments_count,
            runCount: runs_count
        }
    },
    async fetchAllExperiments({page = 1, pageSize = 1000} = {page: 1, pageSize: 1000}) {
        const {data} = await axiosInstance.get(`/epolyscat_django_app/api/experiments/?page=${page}&page_size=${pageSize}`);
        data.results = data.results.map(this.encodeObj);
        return data;
    },
    async fetchExperiment({experimentId}) {
        const {data} = await axiosInstance.get(`/epolyscat_django_app/api/experiments/${experimentId}`);
        return this.encodeObj(data);
    },
    async createExperiment({name, description}) {
        let {data} = await axiosInstance.post("/epolyscat_django_app/api/experiments/", {
            "name": name,
            "description": description
        });

        return this.encodeObj(data);
    },
    async deleteExperiment({experimentId = null} = {}) {
        await axiosInstance.delete(`/epolyscat_django_app/api/experiments/${experimentId}/`);
    }
}

export const PlotService = {
    encodeObj(obj) {
        return {
            xAxis: obj.xaxis,
            yAxes: obj.yaxes,
            flags: obj.flags
        }
    },
 
    async plotSelectedRuns({runIds, expectationValue, xAxis, yAxis, flags}) {
        const {data} = await axiosInstance.post("/epolyscat_django_app/api/plot/", {
            "runs": runIds,
            "plotfiles": plotfiles.map(({ dataProductURI, prefix }) => ({
                "data_product_uri": dataProductURI,
                "prefix": prefix
            })),

            //"plotfile": expectationValue,
            "plot_parameters": {"xaxis": `${xAxis}`, "yaxes": `${yAxis}`, "flags": `${flags}`}
        });

        return {
            "plotImageUrl": data["mime-type"] && data["plot"] ? `data:${data["mime-type"]};base64,${data["plot"]}` : null,
            "output": data["output"],
            "userGuidance": data["user_guidance"]
        }
    },

    async getRunListInputs({runIds}) {
        const {data} = await axiosInstance.post("/epolyscat_django_app/api/list-inputs/", {
            "runs": runIds
        });

        return {
            "output": data.output
        };
    },
    async getRunListInputDifference({runIds}) {
        const {data} = await axiosInstance.post("/epolyscat_django_app/api/diff-inputs/", {
            "runs": runIds
        });

        return {
            "output": data.output
        };
    },
    async getPlotables({runIds}) {
        const {data} = await axiosInstance.post("/epolyscat_django_app/api/plotables/", {
            "runs": runIds
        })

        return data.filenames;
    },
    async getPlotParameters() {
        const { data } = await axiosInstance.get(appBaseUrl + "plot-parameters/");

        return data.map(this.encodeObj);
    },
/*
    async getPlotParameters() {
        let {data} = await axiosInstance.get("/epolyscat_django_app/api/plot-parameters/");
        data = data.map(({xaxis, yaxes, flags}) => {
            return {text: `x=${xaxis} y=${yaxes} ${flags}`, value: {xAxis: xaxis, yAxis: yaxes, flags: flags}};
        });

        return data;
    },
*/
    async getViewables({runId}) {
        const {data} = await axiosInstance.get(`/epolyscat_django_app/api/runs/${runId}/viewables/`)

        return data;
    },
    async getInputFiles({runId}) {
        const {data} = await axiosInstance.get(`/epolyscat_django_app/api/runs/${runId}/input-files/`)

        return data;
    },
}
export const DirectoryService = {
    async fetchDirectories(path = "") {
        const { data } = await axiosInstance.get(`/api/user-storage/${path}`);

        return data;
    },
}
export const SettingsService = {
    async all() {
        const {data} = await axiosInstance.get("/epolyscat_django_app/api/settings/");
        return data;
    },
};

export const RunService = {
    submitAllowedStatuses: ["FAILED", "Unsubmitted"],
    encodeObj(obj) {
        const displayStatus = (
            obj.job_status == "NO STATUS" ||
            obj.job_status == "COMPLETED"
        ) ? obj.status : obj.job_status;
        return {
            runId: obj.id,
            name: obj.name,
            description: obj.description,
            owner: obj.owner,
            id: obj.id,
            airavataProjectId: obj.airavata_project_id,
            jobId: obj.job_id,
            number: obj.number,
            root: obj.root,
            experimentId: obj.experiment,
            created: new Date(obj.created).toLocaleString(),
            updated: new Date(obj.updated).toLocaleString(),
            //created: obj.created,
            //updated: obj.updated,
            deleted: obj.deleted,
            isEmailNotificationOn: obj.is_email_notification_on,
            inpcDownloadUrl: obj.inpc_download_url,

            groupResourceProfileId: obj.group_resource_profile_id,
            computeResourceId: obj.compute_resource_id,
            queueName: obj.queue_name,
            coreCount: obj.core_count,
            nodeCount: obj.node_count,
            wallTimeLimit: obj.walltime_limit,
            totalPhysicalMemory: obj.total_physical_memory,

            status: obj.status,
            jobStatus: obj.job_status,
            resource: obj.resource,
            resourceShort: obj.resource_short,
            executions: obj.executions,
            inputTable: obj.input_table,
            canResubmit: obj.can_resubmit,
            canSubmit: RunService.submitAllowedStatuses.indexOf(obj.status) >= 0,
            inputs: obj.inputs.map(input => ({
                type: input.type,
                name: input.name,
                value: input.value,
                files: input.files.map(file => ({
                    name: file.name,
                    dataProductURI: file.data_product_uri,
                    downloadURL: `sdk/download/?data-product-uri=${encodeURIComponent(file.data_product_uri)}`
                }))
            })),
            viewIds: obj.views,
            displayStatus,
            isTutorial: obj.is_tutorial
        };
      },
    
    async fetchViewableContent({runId, filename, inpcDownloadUrl}) {
        let res;
        if (filename === "inpc" && inpcDownloadUrl) {
            res = await axios.get(inpcDownloadUrl);
        } else {
            res = await axiosInstance.get(`/epolyscat_django_app/api/runs/${runId}/viewables/${filename}/`);
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

        const {data} = await axiosInstance.get("/epolyscat_django_app/api/runs/" + queryString);
        data.results = data.results.map(this.encodeObj);

        return data;
    },
    async fetchRuns() {
        const {data} = await axiosInstance.get(appBaseUrl + 'runs/');

        data.results = data.results.map(this.encodeObj);
        return data;
    },

    async fetchRun({runId = null} = {}) {
        const {data} = await axiosInstance.get(`/epolyscat_django_app/api/runs/${runId}`);
        const _run = this.encodeObj(data);

        if (!_run.inputTable) {
            const _clone = await this.cloneRun({runId});
            _run.inputTable = _clone.inputTable;
        }

        return _run;
    },
    //async cloneRun({runId = null} = {}) {
    //    const {data} = await axiosInstance.post(`/epolyscat_django_app/api/runs/${runId}/new/`);

    //    return this.encodeObj(data);
    //},
    //async deleteRun({runId = null} = {}) {
    //    await axiosInstance.delete(`/epolyscat_django_app/api/runs/${runId}/`);
    //},
    //async createRun({root, experimentId, directedit, inputTable, groupResourceProfileId, computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory, viewIds, description}, submit = false) {
    //async createRun({ root, experimentId, directedit, inputTable, groupResourceProfileId, computeResourceId, queueName, coreCount, nodeCount, wallTimeLimit, totalPhysicalMemory, viewIds, description

    async createRun({ 
        name, inputs, groupResourceProfileId, computeResourceId, coreCount, nodeCount, wallTimeLimit, queueName, totalPhysicalMemory, viewIds, description
    }) {
        let data = {
            "name": name,
            "description": description,
            "inputs_data": inputs, // prevent collision with "inputs" field on the Run serializer       
            //"root": root,
            //"experiment": experimentId,
            //"directedit": directedit,
            //"input_table": inputTable,
            "group_resource_profile_id": groupResourceProfileId,
            "compute_resource_id": computeResourceId,
            "queue_name": queueName,
            "core_count": coreCount,
            "node_count": nodeCount,
            "walltime_limit": wallTimeLimit,
            "total_physical_memory": totalPhysicalMemory,
            viewIds
        };

        //const runCreateRes = await axiosInstance.post("/epolyscat_django_app/api/runs/", data);
        //data = {...data, ...runCreateRes.data};
        const { data: createdRunData} = await axiosInstance.post(appBaseUrl + 'runs/', data);

        return this.encodeObj(createdRunData);

        //if (submit) {
        //    const runSubmitRes = await axiosInstance.post(`/epolyscat_django_app/api/runs/${data.id}/submit/`, data);
        //    data = {...data, ...runSubmitRes.data};
        //}

        //return this.encodeObj(data);
    },

    async updateRun({
        name, runId, inputs, groupResourceProfileId, computeResourceId,
        coreCount, nodeCount, wallTimeLimit, queueName, totalPhysicalMemory, description
    }) {
        let data = {
            name,
            "description": description,
            "inputs_data": inputs, // prevent collision with "inputs" field on the Run serializer
            "group_resource_profile_id": groupResourceProfileId,
            "compute_resource_id": computeResourceId,
            "queue_name": queueName,
            "core_count": coreCount,
            "node_count": nodeCount,
            "walltime_limit": wallTimeLimit,
            "total_physical_memory": totalPhysicalMemory
        };

        const { data: updatedRunData} = await axiosInstance.patch(appBaseUrl + `runs/${runId}/`, data);

        return this.encodeObj(updatedRunData);
    },
/*
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

        const runSubmitRes = await axiosInstance.post(`/epolyscat_django_app/api/runs/${runId}/submit/`, data);
        data = {...data, ...runSubmitRes.data};

        return this.encodeObj(data);
    },
    async updateRun({ 
        name, runId, inputs, groupResourceProfileId, computeResourceId,
        coreCount, nodeCount, wallTimeLimit, queueName, totalPhysicalMemory, description
    }) {
        let data = { 
            name,
            "description": description,
            "inputs_data": inputs, // prevent collision with "inputs" field on the Run serializer
            "group_resource_profile_id": groupResourceProfileId,
            "compute_resource_id": computeResourceId,
            "queue_name": queueName,
            "core_count": coreCount,
            "node_count": nodeCount,
            "walltime_limit": wallTimeLimit,
            "total_physical_memory": totalPhysicalMemory
        };
    
        const { data: updatedRunData} = await axiosInstance.patch(appBaseUrl + `runs/${runId}/`, data);
            
        return this.encodeObj(updatedRunData);
    },      
*/
    async submitRun({ runId }) {
        const { data } = await axiosInstance.patch(appBaseUrl + `runs/${runId}/submit/`, {});
        return this.encodeObj(data);
    },
    async deleteRun({ runId, deleteAssociated }) {
        await axiosInstance.delete(appBaseUrl + `runs/${runId}/?deleteAssociated=${deleteAssociated}`);
    },
    async cloneRun({ runId }) {
        const { data } = await axiosInstance.post(appBaseUrl + `runs/${runId}/clone/`, {});

        return this.encodeObj(data);
    },
    async fetchStatus({runId}) {
        const { data } = await axiosInstance.get(appBaseUrl + `runs/${runId}/status/`, {})

        return data;
    },
    async changeNotificationSettings({ runId, isEmailNotificationOn }) {
        const { data } = await axiosInstance.patch(appBaseUrl + `runs/${runId}/change_notification_settings/`,
            { is_email_notification_on: isEmailNotificationOn }
        );

        // console.log({ is_email_notification_on: isEmailNotificationOn }, data);

        return data;
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

        const runResubmitRes = await axiosInstance.post(`/epolyscat_django_app/api/runs/${runId}/resubmit/`, data);

        data = {...data, ...runResubmitRes.data};

        return this.encodeObj(data);
    },
}

export const ViewService = {
    //readonlyViewTypes: ["unsubmitted", "tutorial"],

    encodeObj(obj) {
        return {
            //viewId: obj.id,
            id: obj.id,
            name: obj.name,
            //owner: obj.owner,
            created: obj.created,
            updated: obj.updated,

            //created: new Date(obj.created).toLocaleString(),
            //updated: new Date(obj.updated).toLocaleString(),
            //deleted: obj.deleted,
            //type: obj.type,
            //activeRunCount: obj.active_run_count,
            runCount: obj.run_count,
            runs: obj.runs.map(RunService.encodeObj)
            //readonly: ViewService.readonlyViewTypes.indexOf(obj.type) >= 0
        };
    },
    async fetchAllViews({page = 1, pageSize = 1000, tutorials = false} = {page: 1, pageSize: 1000, tutorials: false}) {
        let url, data, res;
        if (tutorials) {
            url = "/epolyscat_django_app/api/views/tutorials/";
            res = await axiosInstance.get(url);
            data = {"count": 1, "next": null, "previous": null, "results": [res.data].map(this.encodeObj)}
        } else {
            url = `/epolyscat_django_app/api/views/?page=${page}&page_size=${pageSize}`;
            res = await axiosInstance.get(url);
            data = res.data;
            data.results = data.results.map(this.encodeObj);
        }

        return data;
    },
    async fetchViews() {
        const {data: views} = await axiosInstance.get(appBaseUrl + 'views/');
    
        // return [...views, await this.fetchTutorialsView()].map(this.encodeObj);
        return views.map(this.encodeObj);
    },  
    async createView({
        name, runIds
    }) {
        let data = {
            "name": name,
            "runIds": runIds
        };

        // console.log("Data to send: ", data);

        const { data: view } = await axiosInstance.post(appBaseUrl + 'views/', data);

        return this.encodeObj(view);
    },
    async updateView({
        name, runIds, viewId, overide
    }) {
        let data = {
            "name": name,
            "runIds": runIds,
            "overide": overide
        };

        const { data: view } = await axiosInstance.patch(appBaseUrl + `views/${viewId}/`, data);

        return this.encodeObj(view);
    },
 
    async fetchTutorialsView() {
        const { data: runs } =  await axiosInstance.get(appBaseUrl + `runs/tutorial_runs`);
        
        return this.encodeObj({
            id: -1,
            name: "Tutorials",
            created: (new Date("Nov 3, 2023 10:00:00")).toISOString(),
            updated: (new Date("Nov 3, 2023 10:00:00")).toISOString(),
            run_count: runs.length,
            runs
        });
    },
    async fetchView({viewId}) {
        if (viewId == -1){
           const { data: runs } =  await axiosInstance.get(appBaseUrl + `runs/tutorial_runs`);
           return this.encodeObj({
            id: -1,
            name: "Tutorials",
            created: (new Date("Nov 3, 2023 10:00:00")).toISOString(),
            updated: (new Date("Nov 3, 2023 10:00:00")).toISOString(),
            run_count: runs.length,
            runs
           });
        }
        else{   
           const {data} = await axiosInstance.get(`/epolyscat_django_app/api/views/${viewId}`);
           return this.encodeObj(data);
        }
    },

/*
    async createView({name, runIds}) {
        let {data} = await axiosInstance.post("/epolyscat_django_app/api/views/", {
            "name": name
        });
        const {viewId} = this.encodeObj(data);
        await axiosInstance.put(`/epolyscat_django_app/api/views/${viewId}/add-runs/`, {
            "runs": runIds
        });

        return this.encodeObj({id: viewId, name, runs: runIds});
    },
    async updateView({viewId, name, runIds}) {
        await axiosInstance.put(`/epolyscat_django_app/api/views/${viewId}/`, {
            "name": name
        });
        await axiosInstance.put(`/epolyscat_django_app/api/views/${viewId}/add-runs/`, {
            "runs": runIds
        });

        return this.encodeObj({id: viewId, name, runs: runIds});
    },
*/



    async removeRuns({viewId, runIds}) {
        await axiosInstance.put(`/epolyscat_django_app/api/views/${viewId}/remove-runs/`, {
            "runs": runIds
        });
    },
    async addRuns({viewId, runIds}) {
        await axiosInstance.put(`/epolyscat_django_app/api/views/${viewId}/add-runs/`, {
            "runs": runIds
        });
    },
    async deleteView({viewId = null} = {}) {
        await axiosInstance.delete(`/epolyscat_django_app/api/views/${viewId}/`);
    }
}
