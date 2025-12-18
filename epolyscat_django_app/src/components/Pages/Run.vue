<template>
  <LoadingOverlay class="new-run-page" style="position: relative" :name="loadingKey">
    <!-- div style="position: absolute; top: -1000px;">
            <h3 ref="nameHelper" style="width: fit-content; font-size: 1.75rem"  role="textbox"></h3>
    </div -->
    <!-- 
    <div class="w-100 h-100 bg-light p-2">
     <multipane class="w-100 h-100" layout="vertical">
      <div class="h-100 bg-white overflow-auto p-3" style="flex-grow: 1">
        <div class="w-100">
          <b-breadcrumb>
            <template v-if="experimentId">
              <b-breadcrumb-item to="/experiments">Experiments</b-breadcrumb-item>
              <b-breadcrumb-item v-if="experiment" :to="`/runs?experimentId=${experiment.experimentId}`">
                {{ experiment.name }}
              </b-breadcrumb-item>
            </template>
            <template v-else-if="viewId">
              <b-breadcrumb-item to="/Views">Views</b-breadcrumb-item>
              <b-breadcrumb-item v-if="view" :to="`/runs?viewId=${view.viewId}`">
                {{ view.name }}
              </b-breadcrumb-item>
            </template>
            <b-breadcrumb-item v-else to="/runs">Runs</b-breadcrumb-item>
            <template v-if="runId">
              <b-breadcrumb-item :to="runLink" v-if="run">
                {{ run.name }}
              </b-breadcrumb-item>
            </template>
          </b-breadcrumb>

          <Errors :errors="errors"/>

        </div>
      </div>
    </div>
     -->
    <div style="position: absolute; top: -1000px;">
            <h3 ref="nameHelper" style="width: fit-content; font-size: 1.75rem"  role="textbox"></h3>
    </div>
    <div class="d-flex flex-row align-items-center mb-2 flex-grow-1"></div>

    <div class="w-100 d-flex flex-row justify-content-between" style="padding: 20px 30px 0; flex-wrap: wrap;">
        <div class="d-flex flex-row align-items-center mb-2 flex-grow-1">

                <b-form-input
                    v-if="runViewType != RUN_VIEW_TYPE.viewing" v-model="run.name"
                    style="font-size: 1.75rem" :style="{ width: nameWidth() }"
                />
                <h3 style="margin: 0;" v-else>{{ run.name }}</h3>
                <b-button
                    variant="link" size="sm" @click="moreInfoExpanded = !moreInfoExpanded"
                    :aria-expanded="moreInfoExpanded ? 'true' : 'false'"
                    aria-controls="moreInfo" class="ml-2"
                >
                    <b-icon v-if="moreInfoExpanded" icon="chevron-down"/>
                    <b-icon v-else icon="chevron-right"/>
                </b-button>
        </div>
        <span></span>
        <div>
                <b-button
                    v-if="runViewType != RUN_VIEW_TYPE.creating" variant="outline-primary"
                    @click="cloneRun"
                >Clone</b-button>
                <b-button
                    v-if="run.status.toLocaleUpperCase() == 'COMPLETED'" variant="outline-primary"
                    @click="newRunWithOutputs"
                >Create new run with outputs</b-button>
                <b-button
                    v-if="runViewType != RUN_VIEW_TYPE.creating && !viewOnly" variant="outline-danger"
                    @click="deleteRun"
                >Delete</b-button>
        </div>
    </div>
    <b-collapse v-model="moreInfoExpanded" id="moreInfo" class="m-0 ml-3">
     <div class="d-flex flex-column ml-4 px-3 py-2" style="border-left: 1px solid black; "> 

      <span v-if="runViewType == RUN_VIEW_TYPE.viewing" class="mb-2"><b>Run Type:</b> {{ runType }}</span>
      <div v-if="run.jobId != null" class="mb-2">
         <b>Job ID:</b> {{ run.jobId }}
      </div>
      <b>Description:</b>
               <b-form-textarea
                    id="run-description" v-model="run.description"
                    placeholder="Enter a description for the run"
                    rows="3" style="width: min(600px, 100%)"
                    :readonly="runViewType == RUN_VIEW_TYPE.viewing"
      ></b-form-textarea>
     </div>
    </b-collapse>
    <!-- div class="separator" style="width: calc(100% - 20px)" -->

        <!-- div class="d-flex flex-column" style="margin: 0 40px">
            <span style="width: 75%">
                <b>Descriptions</b> &nbsp;
                <b-button
                    variant="link" @click="showDescriptions = !showDescriptions"
                    style="vertical-align: unset; outline: 0; box-shadow: none;"
                >
                    <span v-if="!showDescriptions">(Show)</span>
                    <span v-else>(Hide)</span>
                </b-button>
            </span>
            <b-collapse visible id="descriptions" v-model="showDescriptions">
                <div style="border-top: 1px solid black; margin-bottom: 10px;">
                </div>
                <div v-for="description in pathDescriptions" :key="description">
                     <b>{{ description.name }}:</b> {{ description.text }}
                </div>
            </b-collapse>
        </div -->

        <div class="separator" style="width: calc(100% - 20px)">

            <b class="ml-2">Experiment status:</b>
            <Badge :status="run.status" class="mx-3" ref="badge"></Badge>
            <b-spinner
                v-if="expStatusSpinner"
                label="Spinning" variant="info" class="mr-3" small
            />
            <b class="ml-2">Job status:</b>
            <Badge :status="run.jobStatus" class="mx-3" ref="badge"></Badge>
            <b-spinner
                v-if="jobStatusSpinner"
                label="Spinning" variant="info" class="mr-3" small
            />
        </div >
        <PathSelecter :viewing="runViewType == RUN_VIEW_TYPE.viewing" ref="pathSelecter" @clickedPath="goToFileOptions"/>
        <div class="d-flex flex-column" style="margin: 0 40px">
            <span style="width: 75%">
                <b>Descriptions</b> &nbsp;
                <b-button
                    variant="link" @click="showDescriptions = !showDescriptions"
                    style="vertical-align: unset; outline: 0; box-shadow: none;"
                >
                    <span v-if="!showDescriptions">(Show)</span>
                    <span v-else>(Hide)</span>
                </b-button>
            </span>
            <b-collapse visible id="descriptions" v-model="showDescriptions">
                <div style="border-top: 1px solid black; margin-bottom: 10px;">
                </div>
                <div v-for="description in pathDescriptions" :key="description">
                    <b>{{ description.name }}:</b> {{ description.text }}
                </div>
            </b-collapse>
        </div>
        <FileOptions :viewType="viewType" :runViewType="runViewType" :runId="runId" ref="fileOptions"/>
        <h3>Parameters</h3>
        <div style="padding: 0 20px; width: 100%">
            <TableView v-show="parameters.length > 0" :key="paramTableKey" :viewing="runViewType == RUN_VIEW_TYPE.viewing" :selectedTableObject="parametersTableObject" ref="parametersTable"/>
        </div>
        <p v-show="parameters.length == 0">No parameters</p>
        <RunResource :run="run" :viewing="runViewType == RUN_VIEW_TYPE.viewing" @updateResources="updateResources"/>
         <div class="w-100 d-flex flex-row justify-content-between createRunButtons" v-if="runViewType != RUN_VIEW_TYPE.viewing">
            <!-- Wrapper to allow for tooltips to show when disabled -->
            <span id="saveButton" style="width: fit-content">
                <b-button variant="primary" :disabled="saveIssues.length > 0" @click="saveRun()">Save</b-button>
            </span>
            <b-tooltip target="saveButton" placement="top" v-if="saveIssues.length > 0">
                <div v-html="saveIssues.join('<br>')" />
            </b-tooltip>

            <!-- Wrapper to allow for tooltips to show when disabled -->
            <span id="submitButton" style="width: fit-content">
                <b-button variant="primary" :disabled="submitIssues.length > 0" @click="submitRun">Save & Submit</b-button>
            </span>
            <b-tooltip target="submitButton" placement="top" v-if="submitIssues.length > 0">
                <div v-html="submitIssues.join('<br>')" />
            </b-tooltip>
         </div>
  </LoadingOverlay>
</template>

<script>
//import RunActions from "@/components/block/RunActions";
//import {PlotService} from "@/service/epolyscat-service";
//import Errors from "@/components/block/Errors";
//import {Multipane, MultipaneResizer} from "vue-multipane";
//import PostProcessing from "@/components/block/PostProcessing";
//import RunViewableEditor from "@/components/block/RunViewableEditor";

import PathSelecter from "@/components/blocks/PathSelecter";
import UserStorage from "@/components/overlay/UserStorage";
import RunResource from "@/components/blocks/RunResource";
import TableView from "@/components/blocks/TableView";
import FileOptions from "../block/FileOptions";
import Badge from "../block/Badge.vue";
import store from "@/store";
import router from '@/router';
import { eventBus } from "@/event-bus";
import LoadingOverlay from "@/components/overlay/LoadingOverlay";
import { descriptions, tableObjects } from "@/fileData";


export default {
  store,
  router,
  created() {
        this.VIEW_TYPE = {
            none: "none",
            upload: "upload",
            text: "text",
            table: "table"
        };

        this.RUN_VIEW_TYPE = {
            creating: "creating",
            editing: "editing",
            viewing: "viewing"
        }
  },
  //name: 'CreateRun',
  //components: {RunResource, RunViewableEditor, PostProcessing, Errors, RunActions, Multipane, MultipaneResizer},
  //store: store,
  data() {
    let moreInfoExpanded = this.$store.getters["settings/getPreference"]("moreInfoExpanded");
    let showDescriptions = this.$store.getters["settings/getPreference"]("showDescriptions");
    if (moreInfoExpanded == null) moreInfoExpanded = false;
    if (showDescriptions == null) showDescriptions = false;

    return {
      //plotables: null,
      //viewables: null,
      //filename: "inpc",
      //inputFiles: null,
      //processingPlotables: null,
      //processingViewables: null,
      //processingInputFiles: null,
      //processing: false,
      //selectedPlotable: null,
      runId: -1, // -1 means that the run has yet to be created and any other value is a real runId
      //errors: [],
      userLeftPage: true,
      moreInfoExpanded,
      showDescriptions,
      spinnerStates: [
             "CREATED", "SUBMITTED", "VALIDATED", "QUEUED", "SCHEDULED",
             "LAUNCHED", "EXECUTING", "ACTIVE", "CANCELING"
      ],
      paramTableKey: 0,
      viewIds: [],
      run: {
          name: `EPolyScat Run on ${(new Date(Date.now())).toLocaleString()}`,
          id: -1,
          status: "UNSUBMITTED",
          jobStatus: "UNSUBMITTED",
          displayStatus: "UNSUBMITTED",
          isEmailNotificationOn: false,
          description: "",
          groupResourceProfileId: null,
          computeResourceId: null,
          queueName: null,
          coreCount: null,
          nodeCount: null,
          wallTimeLimit: null,
          totalPhysicalMemory: null,
          executions: [],
          viewIds: []
      },
      //runRefreshTimeout: null
    };
  },
  components: { 
        PathSelecter, UserStorage, RunResource, TableView, 
        LoadingOverlay, FileOptions, Badge
  },  

  computed: {
    runViewType() {
        return (this.runId == -1) ? this.RUN_VIEW_TYPE.creating
            : (this.run.id == -1) ? this.RUN_VIEW_TYPE.viewing
            : ((this.run.executions == undefined || this.run.executions.length == 0) && !this.run.isTutorial) ? this.RUN_VIEW_TYPE.editing
            : this.RUN_VIEW_TYPE.viewing;
    },
    viewOnly() {
            return this.run.isTutorial;
    },
    viewType() {
        const { selected, isFile } = this.$store.getters["input/getSelectedInfo"]();

        if (selected == undefined)
             return this.VIEW_TYPE.none;

        if (selected == "Upload")
             return this.VIEW_TYPE.upload;

        if (isFile)
             return this.VIEW_TYPE.text;

        if (selected == "Parameters")
             return this.VIEW_TYPE.table;

        return this.VIEW_TYPE.none;
    },
    saveIssues() {
            return [];
    },
    submitIssues() {
            const { inputFiles, parameters } = this.$store.getters["input/getInputs"]();

            let issues = Object.entries(inputFiles)
                    .filter(([_, inputFile]) => !inputFile.isValid())
                    .map(([filename, _]) => `"${filename}": still needs files to be uploaded`)
                .concat(parameters
                    .filter(parameter => !parameter.issues(parameter.value).length == 0)
                    .map(parameter => `"${parameter.name}": "${parameter.value}" is not valid`)
                );

            if (this.run.groupResourceProfileId == null) issues.push("The allocation still needs to be set");
            if (this.run.computeResourceId == null) issues.push("The compute resource still needs to be set");
            if (this.run.queueName == null) issues.push("The queue still needs to be set");
            if (this.run.coreCount == null) issues.push("The number of cores still needs to be set");
            if (this.run.nodeCount == null) issues.push("The number of nodes still needs to be set");
            if (this.run.wallTimeLimit == null) issues.push("The wall time limit still needs to be set");
            if (this.run.totalPhysicalMemory == null) issues.push("The total physics memory still needs to be set");

            return issues
    },
    pathDescriptions() {
            const path = this.$store.getters["input/getPath"]
                .slice((this.runViewType == this.RUN_VIEW_TYPE.viewing) ? 2 : 0)
                .filter(item => item in descriptions);

            return path
                .filter((item, i) => path.indexOf(item) == i)
                .map(item => ({
                    name: item,
                    text: descriptions[item]
                }));
    },
    loadingKey() {
            return `page-${this.run.name}`;
    },
    runType() {
            if (this.run.inputs == null || this.run.inputs.length == 0)
                return "-----";

            let type1 = this.run.inputs.find(input => input.name == "Calculation_Type");
            let type2 = this.run.inputs.find(input =>
                ["EPOLYSCAT_Application_Module", "Application_Utility", "Application_Workflow"].indexOf(input.name) != -1
    );

            return `${type1.value} > ${type2.value}`;
    },
    parameters() {
            return this.$store.getters["input/getInputs"]().parameters;
    },
    parametersTableObject() {
            return tableObjects["Parameters"];
    },
    expStatusSpinner() {
            const executingStates = [
                "CREATED", "SUBMITTED", "VALIDATED", "QUEUED", "SCHEDULED",
                "LAUNCHED", "EXECUTING", "ACTIVE", "CANCELING"
            ];

            return executingStates.indexOf(this.run.status.toLocaleUpperCase()) != -1;
    },
    jobStatusSpinner() {
            const executingStates = [
                "CREATED", "SUBMITTED", "VALIDATED", "QUEUED", "SCHEDULED",
                "LAUNCHED", "EXECUTING", "ACTIVE", "CANCELING"
            ];

            return executingStates.indexOf(this.run.jobStatus.toLocaleUpperCase()) != -1;
    },


/*

    filenames() {
      if (this.viewables) {
        return this.viewables.map(({filename}) => filename);
      } else {
        return null;
      }
    },
//    runLink() {
//      let _link = "/runs/";
//      if (this.runId) {
//        _link += `${this.runId}?`;
//      }
//
//      if (this.experimentId) {
//        _link += `experimentId=${this.experimentId}&`;
//      }
//
//      if (this.viewId) {
//        _link += `viewId=${this.viewId}&`;
//      }
//
//      return _link;
//    },
    experimentId() {
      return this.$route.query.experimentId;
    },
    viewId() {
      return this.$route.query.viewId;
    },
    runId() {
      return this.$route.params.runId;
    },
    run() {
      return this.$store.getters["run/getRun"]({
        runId: this.runId
      });
    },
    experiment() {
      return this.$store.getters["experiment/getExperiment"]({
        experimentId: this.experimentId
      });
    },
    view() {
      return this.$store.getters["view/getView"]({
        viewId: this.viewId
      });
    },
    epolyscatApplicationModuleId() {
      return this.$store.getters["settings/epolyscatApplicationModuleId"];
    },
    expStatusSpinner() {
        const executingStates = [
            "CREATED", "SUBMITTED", "VALIDATED", "QUEUED", "SCHEDULED",
            "LAUNCHED", "EXECUTING", "ACTIVE", "CANCELING"
        ];

        return executingStates.indexOf(this.run.status.toLocaleUpperCase()) != -1;
    },
    jobStatusSpinner() {
        const executingStates = [
                "CREATED", "SUBMITTED", "VALIDATED", "QUEUED", "SCHEDULED",
                "LAUNCHED", "EXECUTING", "ACTIVE", "CANCELING"
        ];

        return executingStates.indexOf(this.run.jobStatus.toLocaleUpperCase()) != -1;
    },
*/
  },
  methods: {
      async fetchRun() {
            if (this.runId != -1) {
                this.$store.commit("loading/START", { key: this.loadingKey, message: "Fetching Run" });

                try {
                    await this.$store.dispatch("run/fetchRun", { runId: this.runId });
                    await this.$store.dispatch("run/loadInputs", { runId: this.runId });
                    this.run = await this.$store.dispatch("run/tryLoadOutputs", { runId: this.runId });

                    if (!this.userLeftPage)
                        this.$refs.pathSelecter.selectDefaultPathFrom(2);

                    // await this.$store.dispatch("input/fetchOutputs", { runId: this.runId });
                    this.$store.commit("loading/STOP", { key: this.loadingKey, message: "Fetching Run" });
                } catch (error) {
                    eventBus.$emit("error", { name: `Error while trying to fetch run with id: ${this.runId}`, error });
                }
            }

            // console.log("the run:", {...this.run})
      },
      async fetchURLData() {
            if ("viewId" in this.$route.query) {
                this.viewIds = [this.$route.query.viewId];
            } else if ("withOutputsFrom" in this.$route.query) {
                const outputsRunId = parseInt(this.$route.query.withOutputsFrom);

                this.$store.commit("loading/START", { key: this.loadingKey, message: "Fetching Outputs" });

                try {
                    const run = await this.$store.dispatch("run/fetchRun", { runId: outputsRunId });
                    await this.$store.dispatch("input/loadOutputsIntoInputs", { runId: outputsRunId });
                    this.run.name = `Continuation of ${run.name}`;
                    this.viewIds = (run.isTutorial) ? [] : run.viewIds;

                    this.$store.commit("loading/STOP", { key: this.loadingKey, message: "Fetching Outputs" });
                } catch (error) {
                    eventBus.$emit("error", {
                        name: `Error while trying to fetch the outputs from the run with id: ${outputsRunId}`,
                        error
                    });
                }
            } else if ("clonedFrom" in this.$route.query) {
                const clonedRunId = parseInt(this.$route.query.clonedFrom);

                this.$store.commit("loading/START", { key: this.loadingKey, message: "Fetching Cloned Run" });

                try {
                    const run = await this.$store.dispatch("run/fetchRun", { runId: clonedRunId });
                    await this.$store.dispatch("run/loadInputs", { runId: clonedRunId });
                    this.viewIds = (run.isTutorial) ? [] : run.viewIds;
                    const strippedRun = this.$store.getters["run/getStrippedRun"](run);

                    this.run = {
                        ...this.run,
                        ...strippedRun,
                        name: `Clone of ${run.name}`,
                        id: -1,
                    };
                    // console.log("run: ", run, "stripped: ", strippedRun, " this.run: ", this.run);
                    if (!this.userLeftPage)
                        this.$refs.pathSelecter.selectDefaultPathFrom(2);

                    this.$store.commit("loading/STOP", { key: this.loadingKey, message: "Fetching Cloned Run" });
                } catch (error) {
                    eventBus.$emit("error", {
                        name: `Error while trying to clone run with id: ${clonedRunId}`,
                        error
                    });
                }
            }
      },
      async saveRun(setURL = true) {
            this.$store.commit("loading/START", { key: this.loadingKey, message: "Saving Run" });

            if (this.runViewType == "creating") {
                try {
                    this.run = await this.$store.dispatch("run/createRun", {
                        ...this.run, viewIds: this.viewIds
                    });

                    if (!this.userLeftPage && setURL) router.push(`/runs/${this.run.id}`)
                } catch (error) {
                    eventBus.$emit("error", { name: `Error while trying to create the run`, error });
                }

                this.runId = this.run.id;
            } else {
                try {
                    this.run = await this.$store.dispatch("run/updateRun", {
                        ...this.run
                    });
                } catch (error) {
                    eventBus.$emit("error", { name: `Error while trying to update the run`, error });
                }
            }

            this.$store.commit("loading/STOP", { key: this.loadingKey, message: "Saving Run" });
      },
      async submitRun() {
            await this.saveRun(false);

            this.$store.commit("loading/START", { key: this.loadingKey, message: "Submitting Run" });

            try {
                this.run = await this.$store.dispatch("run/submitRun", {
                    runId: this.runId
                });
            } catch (error) {
                eventBus.$emit("error", { name: `Error while trying to submit the run`, error });
            }
            if (!this.userLeftPage)
                this.$refs.pathSelecter.selectDefaultPathFrom(2);
            // await this.$store.dispatch("input/fetchOutputs", { runId: this.runId });

            this.$store.commit("loading/STOP", { key: this.loadingKey, message: "Submitting Run" });

            // await this.fetchRun();

            this.run.status = "EXECUTING";
            this.run.jobStatus = "NO STATUS";

            if (!this.userLeftPage) router.push(`/runs/${this.run.id}`)
      },
      async cloneRun() {
            if (!this.userLeftPage) {
                await this.$router.push(`/runs/new?clonedFrom=${this.run.id}`);
                this.$router.go()
            }
            // this.$store.commit("loading/START", { key: this.loadingKey, message: `Cloning run` });

            // try {
            //     const run = await this.$store.dispatch("run/cloneRun", this.run);

            //     this.$store.commit("loading/STOP", { key: this.loadingKey, message: `Cloning run` });
            //     if (!this.userLeftPage) {
            //         await this.$router.push(`/runs/${run.id}`);
            //         this.$router.go()
            //     }
            // } catch (error) {
            //     eventBus.$emit("error", { name: `Error while trying to clone the run with id: ${this.run.id}`, error });
            // }

      },
      async newRunWithOutputs() {
            if (!this.userLeftPage) {
                await this.$router.push(`/runs/new?withOutputsFrom=${this.run.id}`);
                this.$router.go()
            }
      },
      async deleteRun() {
            const deleteAll = await this.$bvModal.msgBoxConfirm(
                `Are you sure you want to delete: "${this.run.name}"? (This will redirect you away from this page)`, {
                    title: 'Please Confirm',
                    okVariant: 'danger',
                    okTitle: 'DELETE',
                    cancelTitle: 'Cancel',
                    footerClass: 'p-2',
                    hideHeaderClose: false,
                    centered: true
                }
            );

            if (deleteAll) {
                const deleteAssociated = await this.$bvModal.msgBoxConfirm(
                    `Do want to delete all files associated with "${this.run.name}"?`, {
                        title: 'Please Confirm',
                        okVariant: 'danger',
                        okTitle: 'Yes',
                        cancelTitle: 'No',
                        footerClass: 'p-2',
                        hideHeaderClose: false,
                        centered: true
                    }
                )

                this.$store.commit("loading/START", { key: this.loadingKey, message: `Deleting run` });

                try {
                    await this.$store.dispatch("run/deleteRun", { runId: this.run.id, deleteAssociated 
                    });

                    this.$store.commit("loading/STOP", { key: this.loadingKey, message: `Deleting run` });
                    this.$router.push(`/`);
                } catch (error) {
                    eventBus.$emit("error", {
                        name: `Error while trying to delete "${this.run.name}"`,
                        error
                    });
                }
            }
      },
      updateResources(updatedResources) {
            console.log("UPDATING RUN FROM RESOURCES: ", this.run, {...this.run})
            this.run.groupResourceProfileId = updatedResources.groupResourceProfileId;
            this.run.computeResourceId = updatedResources.computeResourceId;
            this.run.nodeCount = updatedResources.nodeCount;
            this.run.queueName = updatedResources.queueName;
            this.run.coreCount = updatedResources.coreCount;
            this.run.wallTimeLimit = updatedResources.wallTimeLimit;
            this.run.totalPhysicalMemory = updatedResources.totalPhysicalMemory;
      },
      nameWidth() {
            if ("nameHelper" in this.$refs) {
                this.$refs.nameHelper.innerHTML = this.run.name.replaceAll(" ", "&nbsp;");
                let width = this.$refs.nameHelper.clientWidth + 24;

                return `min(${width}px, 100%)`
            } else {
                return `0px`
            }
      },
      goToFileOptions() {
            this.$refs.fileOptions.$el.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
        },
        async updateNotificationSettings() {
            try {
                this.run = this.$store.dispatch("run/changeNotificationSettings", {
                    runId: this.runId,
                    isEmailNotificationOn: this.run.isEmailNotificationOn
                });
            } catch (error) {
                eventBus.$emit("error", { name: `Error while trying to update notification settings`, error });
            }
      },



      async refreshData() {
            let runHasBeenSubmitted = (this.run.executions.length > 0) || this.run.status != "UNSUBMITTED";
            let runIsntFinished = ["COMPLETED", "FAILED"]
                .indexOf(this.run.status.toLocaleUpperCase()) == -1;

            if (runHasBeenSubmitted && runIsntFinished) {
                const { pathLabels: pathLabelsBefore } = this.$store.getters["input/getInputs"]();

                try {
                    let run = await this.$store.dispatch("run/fetchRun", { runId: this.runId });
                    await this.$store.dispatch("run/loadInputs", { runId: this.runId });
                    this.run = await this.$store.dispatch("run/tryLoadOutputs", { runId: this.runId });

                    const { pathLabels: pathLabelsAfter } = this.$store.getters["input/getInputs"]();

                    if (!this.userLeftPage && pathLabelsBefore != pathLabelsAfter)
                        this.$refs.pathSelecter.selectDefaultPathFrom(2);
                } catch (error) {
                    eventBus.$emit("error", { name: `Error while trying to fetch run with id: ${this.runId}`, error });
                }
            }
      },
    

/*
    viewableLink({filename}) {
      let _link = `/runs/${this.runId}/viewables/${filename}?`;

      if (this.experimentId) {
        _link += `experimentId=${this.experimentId}&`;
      }

      if (this.viewId) {
        _link += `viewId=${this.viewId}&`;
      }

      return _link;
    },
    onDelete() {
      let redirectUrl = "/runs?";
      if (this.experimentId) {
        redirectUrl += `experimentId=${this.experimentId}&`;
      }

      if (this.viewId) {
        redirectUrl += `viewId=${this.viewId}&`;
      }

      this.$router.history.push(redirectUrl);
    },
    refreshData() {
      this.errors = [];

      if (this.runId) {
        this.$store.dispatch("run/fetchRun", {runId: this.runId});

        PlotService.getPlotables({runIds: [this.runId]}).then(plotables => {
          this.plotables = plotables;
          this.processingPlotables = false;
        }).catch(e => {
          this.processingPlotables = false;
          this.errors.push({
            variant: "danger",
            title: "Network Error",
            description: "Unknown error when fetching plotables.",
            source: e
          });
        });

        PlotService.getViewables({runId: this.runId}).then(viewables => {
          this.viewables = viewables;
          this.processingViewables = false;
        }).catch(e => {
          this.processingViewables = false;
          this.errors.push({
            variant: "danger",
            title: "Network Error",
            description: "Unknown error when fetching viewables.",
            source: e
          });
        });

        PlotService.getInputFiles({runId: this.runId}).then(inputFiles => {
          this.inputFiles = inputFiles;
          this.processingInputFiles = false;
        }).catch(e => {
          this.processingInputFiles = false;
          this.errors.push({
            variant: "danger",
            title: "Network Error",
            description: "Unknown error when fetching input files.",
            source: e
          });
        });
      }

      if (this.experimentId) {
        this.$store.dispatch("experiment/fetchExperiment", {experimentId: this.experimentId});
      }

      if (this.viewId) {
        this.$store.dispatch("view/fetchView", {viewId: this.viewId});
      }
    },
    setPlotableSelected(plotable) {
      this.selectedPlotable = plotable;
    },
    togglePlotableSelection(plotable) {
      if (this.selectedPlotable !== plotable) {
        this.setPlotableSelected(plotable);
      } else {
        // Toggle the selection.
        this.selectedPlotable = null;
      }
    }
*/
  },
  watch: {
        moreInfoExpanded() {
            this.$store.commit("settings/setPreference", {
                key: "moreInfoExpanded",
                value: this.moreInfoExpanded
            });
        },
        showDescriptions() {
            this.$store.commit("settings/setPreference", {
                key: "showDescriptions",
                value: this.showDescriptions
            });
        },
        parameters() {
            this.paramTableKey++;
        },
        // run(newRun, oldRun) {
        //     console.log("RUNCHANGED: ", {...oldRun}, "TO: ", {...newRun})
        // }
  },

  mounted() {
    this.errors = [];
    if ("runId" in this.$route.params)
            this.runId = this.$route.params.runId;

    this.fetchRun();
    this.fetchURLData();
/*
    this.plotables = null;
    this.processingPlotables = true;

    this.viewables = null;
    this.processingViewables = true;

    this.inputFiles = null;
    this.processingInputFiles = true;

    this.refreshData();

    this.runRefreshTimeout = setInterval(() => {
      this.refreshData();
    }, 15000);
*/
    this.refreshDataInterval = setInterval(this.refreshData, 15000);
  },
  beforeDestroy() {
    clearInterval(this.runRefreshTimeout);
    // console.log("NO MORE REDIRECTING!!!\n\n\n\n");
    this.userLeftPage = true;
    this.$store.commit("loading/CLEAR", { key: this.loadingKey });
  }
};
</script>

<style scoped>
    .new-run-page {
        display: flex;
        flex-direction: column;
        align-items: center;
        width: 100%;
        height: 100%;
    }

    .new-run-page > * {
        margin: 10px 20px;
        width: -webkit-fill-available;
        width:-moz-available;
    }

    .inputs > * {
        margin: 0 10px;
    }

    /* buttons wrappers */
    .createRunButtons > * {
        display: flex;
        flex-grow: 1;
        margin: 20px;
    }

    /* buttons */
    .createRunButtons > * > * {
        flex-grow: 1;
    }
    .separator {
        display: flex;
        align-items: center;
        text-align: center;
    }

    .separator::before,
    .separator::after {
        content: '';
        flex: 1;
        border-bottom: 1px solid black;
    }

    .separator:not(:empty)::before {
        margin-right: .25em;
    }

    .separator:not(:empty)::after {
        margin-left: .25em;
    }
/*



.multipane .multipane-resizer {
  margin: 0;
  left: 0;
  position: relative;
}

.multipane .multipane-resizer:before {
  display: block;
  content: "";
  width: 90%;
  height: 40px;
  position: absolute;
  top: 50%;
  left: 5%;
  margin-top: -20px;
  margin-left: -10%;
  background-color: #237abb;
  border-radius: 10px;
}

.multipane .multipane-resizer:hover,
.multipane .multipane-resizer:before {
  border-color: #999;
}
*/
</style>
