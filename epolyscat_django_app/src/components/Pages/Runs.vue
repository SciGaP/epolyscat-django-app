<template>
  <div class="w-100 h-100 bg-light p-2">
    <multipane class="w-100 h-100" layout="vertical">
      <div class="h-100 bg-white d-flex flex-column overflow-auto align-items-center" style="flex-grow: 1">
                <div class="d-flex flex-row justify-content-between" style="padding: 30px 0 20px; flex-wrap: wrap; width: 93%">
                    <h3 style="line-height: 38px; margin: 0;">{{ pageName }}</h3>
                    <span></span>
                    <div class="d-flex flex-row button_row align-items-center">
                        <span v-if="numberOfRunsSelected > 0" style="line-height: 41px;">
                            {{ numberOfRunsSelected }}
                            <span v-if="numberOfRunsSelected == 1">Run</span>
                            <span v-else>Runs</span>
                            selected
                        </span>
                        <b-button v-if="!isTutorials && numberOfRunsSelected == 1" variant="danger" @click="deleteRuns(selected)">Delete Run</b-button>
                        <b-button v-if="!isTutorials && numberOfRunsSelected > 1" variant="danger" @click="deleteRuns(selected)">Delete Runs</b-button>
                        <b-button v-if="!isTutorials && numberOfRunsSelected > 0 && view != null" variant="danger" @click="removeFromView(selected)">Remove from View</b-button>
                        <b-button id="save-selection-button" v-show="numberOfRunsSelected > 0" variant="light">
                            <b-icon icon="download"/> Save selection
                        </b-button>
                        <b-popover
                            target="save-selection-button"
                            ref="saveSelectionPopover"
                            triggers="click "
                            placement="auto"
                        >
                            <template #title>
                                <div class="d-flex flex-row">
                                    <div class="flex-fill" style="line-height: 27px;">Save selection</div>
                                    <div>
                                        <b-button variant="link" size="sm"
                                            v-on:click="$refs.saveSelectionPopover.$emit('close')"
                                        ><b-icon icon="x"/></b-button>
                                    </div>
                                </div>
                            </template>
                            <SaveSelectedRunsToViewForm
                                :runIds="selected.map(run => run.id)"
                                @close="$refs.saveSelectionPopover.$emit('close')"
                            />
                        </b-popover>
                        <div v-show="!isTutorials && numberOfRunsSelected >= 1" class="mx-1" style="width: 0.1px; height: 50%; border: 0.5px solid black" />
                        <router-link v-if="!isTutorials" :to="newRunLink" v-slot="{ href, navigate }">
                            <b-button variant="light" :href="href" @click="navigate">New Run</b-button>
                        </router-link>
                    </div>




        <!-- div class="w-100 p-3">

          <b-breadcrumb>
            <template v-if="experimentId">
              <b-breadcrumb-item to="/experiments">Experiments</b-breadcrumb-item>
              <b-breadcrumb-item v-if="experiment" :to="`/runs/?experimentId=${experimentId}`">
                {{ experiment.name }}
              </b-breadcrumb-item>
            </template>
            <template v-else-if="viewId">
              <b-breadcrumb-item to="/views">Views</b-breadcrumb-item>
              <b-breadcrumb-item v-if="view" :to="`/runs/?viewId=${viewId}`">
                {{ view.name }}
              </b-breadcrumb-item>
            </template>
          </b-breadcrumb>
        </dev -->
          <div class="w-100 d-flex">
            <div style="flex: 1;" class="overflow-auto">
              <div class="d-inline mr-2">
                <div class="overflow-auto" style="font-weight: 400; font-size: 19px;">
                  <template v-if="!!experiment">{{ experiment.name }}</template>
                  <template v-if="!!view">{{ view.name }}</template>
                  <template v-if="!experiment && !view">Runs</template>
                </div>
                <div class="overflow-auto" v-if="!!experiment" style="font-weight: 300; font-size: 14px;">
                  {{ experiment.description }}
                </div>
              </div>
              <b-badge v-if="!!experiment">{{ experiment.status }}</b-badge>
            </div>
            <div class="pl-3" v-if="!view || !view.readonly">
                 <div v-show="!isTutorials && numberOfRunsSelected >= 1" class="mx-1" style="width: 0.1px; height: 50%; border: 0.5px solid black;">
                   <router-link v-if="!isTutorials" :to="newRunLink" v-slot="{ href, navigate }">
                    <b-button variant="light" :href="href" @click="navigate">New Run</b-button>
                   </router-link>
                 </div>
            </div>

              <LoadingOverlay :name="loadingOverlayName" class="d-flex w-100" style="height: calc(100% - 145px)">
                    <ListView
                        :items="runs" :columns="[['name', 'Name'], ['status', 'Status'], ['updated', 'Last Modified'], ['actions', 'Actions']]"
                        :canSelectMultiple="true" @updateSelected="updateSelected" identifier="id"
                        :sorters="[
                            (run1, run2) => run1.name.localeCompare(run2.name),
                            (run1, run2) => sortStatus(run1.displayStatus, run2.displayStatus),
                            (run1, run2) => (new Date(run2.updated)).getTime() - (new Date(run1.updated)).getTime(),
                        ]"
                    >
                        <template v-slot:name="{ item }">
                            <router-link :to="`/runs/${item.id}`" v-slot="{isExactActive, href, navigate}">
                                <b-link variant="link" :class="{active: isExactActive}" :href="href" @click="navigate">
                                    {{ item.name }}
                                </b-link>
                            </router-link>
                        </template>
                        <template v-slot:status="{ item }">
                            <Badge :status="item.displayStatus"/>
                        </template>
                        <template v-slot:updated="{ item }">
                            {{ (new Date(item.updated)).toLocaleString() }}
                        </template>
                        <template v-slot:actions="{ item }">
                            <b-button
                                v-if="!isTutorials && view != null" variant="link" size="sm" @click="removeFromView([item])"
                                v-b-tooltip.hover title="Remove from view"
                            ><b-icon icon="x" /></b-button>
                            <b-button
                                variant="link" size="sm" @click="cloneRun(item)"
                                v-b-tooltip.hover title="Clone"
                            ><b-icon icon="clipboard" /></b-button>
                            <b-button
                                v-if="!isTutorials && item.displayStatus.toUpperCase() == 'UNSUBMITTED'" variant="link" size="sm"
                                @click="submitRun(item)" v-b-tooltip.hover title="Submit"
                            ><b-icon icon="arrow-bar-up"/></b-button>
                            <b-button
                                v-else-if="!isTutorials" :disabled="['COMPLETED', 'FAILED'].indexOf(item.status.toUpperCase()) == -1" variant="link"
                                size="sm" @click="resubmitRun(item)"
                                v-b-tooltip.hover title="Resubmit"
                            ><b-icon icon="arrow90deg-right"/></b-button>
                            <b-button v-if="!isTutorials"
                                variant="link" size="sm" @click="deleteRuns([item])"
                                v-b-tooltip.hover title="Delete"
                            ><b-icon icon="trash" /></b-button>
                        </template>
                    </ListView>
              </LoadingOverlay>

<!--

              <button-overlay :show="processingDeleteSelected">
                <b-button variant="outline-primary" size="sm" v-if="selectedCount > 0"
                          v-b-modal="'modal-confirmation-delete-selected-runs'">
                  <b-icon icon="trash"></b-icon>
                  Destroy
                </b-button>
              </button-overlay>
              <b-modal id="modal-confirmation-delete-selected-runs" title="Delete Confirmation" ok-title="Delete"
                       v-on:ok="deleteAllSelectedRuns">
                <p class="my-4">Are you sure that you want to delete {{ selectedCount }} selected runs? </p>
              </b-modal>
          </div>
        </div>

          <div class="w-100 mt-3">
            <b-form-input size="sm" type="text" placeholder="Search a run"/>
          </div>

        <div class="w-100 p-3 overflow-auto">
          <table-overlay-info :data="runs" :rows="5" :columns="5" empty-label="No runs available.">
            <b-table-simple>
              <b-thead>
                <b-tr>
                  <b-th></b-th>
                  <b-th>Name</b-th>
                  <b-th v-if="!view || view.type !== 'tutorial'">Exp. Status</b-th>
                  <b-th v-if="!view || view.type !== 'tutorial'">Job Status</b-th>
                  <b-th v-if="!view || view.type !== 'tutorial'">Resource</b-th>
                  <b-th>Actions</b-th>
                </b-tr>
              </b-thead>
              <b-tbody>
                <b-tr v-for="run in runs" :key="run.runId">
                  <b-td>
                    <b-form-checkbox v-model="selectedRunIdsMap[run.runId]"/>
                  </b-td>
                  <b-td>
                    <router-link :to="runLink(run)"
                                 v-slot="{ href, route, navigate, isActive,isExactActive }">
                      <b-button variant="link" size="sm" @click="navigate" v-b-tooltip.hover.auto :title="run.name">
                        <div class="overflow-auto" style="max-width: 200px;">{{ run.name }}</div>
                      </b-button>
                    </router-link>
                  </b-td>
                  <b-td v-if="!view || view.type !== 'tutorial'">{{ run.status }}</b-td>
                  <b-td v-if="!view || view.type !== 'tutorial'">{{ run.JobStatus }}</b-td>
                  <b-td v-if="!view || view.type !== 'tutorial'">
                    <div class="overflow-auto" style="max-width: 200px;">{{ run.resource }}</div>
                  </b-td>
                  <b-td style="min-width: 150px; max-width: 150px;">
                    <RunActions :run="run" :experiment="experiment" :view="view" v-on:delete="refreshList"/>
                  </b-td>
                </b-tr>
              </b-tbody>
            </b-table-simple>
          </table-overlay-info>
        </div>

        <b-pagination
            class="w-100 p-3"
            v-if="runsPagination && runsPagination.total > pageSize"
            v-model="page"
            :total-rows="runsPagination.total"
            :per-page="pageSize">
        </b-pagination>

        <div class="w-100 d-flex p-3" style="border-top: 1px solid #c9c9c9;" v-if="selectedCount">
          <div class="flex-fill">{{ selectedCount }} Runs Selected</div>
          <div v-if="!view || view.type !== 'tutorial'">
            <b-button id="save-selection-button" variant="outline-primary" size="sm" :disabled="selectedCount === 0">
              <b-icon icon="download"/>
              Save selection
            </b-button>
            <b-popover
                target="save-selection-button"
                triggers="click focus"
                placement="auto"
                container="my-container"
                ref="saveSelectionPopover"
                title="Save selection"
            >
              <template #title>
                <div class="d-flex flex-row">
                  <div class="flex-fill" style="line-height: 27px;">Save selection</div>
                  <div>
                    <b-button variant="link" size="sm"
                              v-on:click="$refs.saveSelectionPopover.$emit('close')">
                      <b-icon icon="x"/>
                    </b-button>
                  </div>
                </div>

              </template>
              <save-selected-runs-to-view-form
                  :run-ids="selectedRunIds"
                  @close="$refs.saveSelectionPopover.$emit('close')"/>
            </b-popover>
          </div>
        </div>

-->

      </div>
      <multipane-resizer v-if="isCompareEnabled"/>
      <div class="h-100 m-1" v-if="isCompareEnabled" style="overflow-y: scroll;">
           <ComparePlots :selectedRuns="selected"/>
      </div>

    </multipane>
  </div>
</template>

<script>
import store from "@/store";
//import TableOverlayInfo from "@/components/overlay/table-overlay-info";
import {RunService, ViewService} from "@/service/epolyscat-service";
//import ButtonOverlay from "@/components/overlay/button-overlay";
import {Multipane, MultipaneResizer} from 'vue-multipane';
//import PostProcessing from "@/components/block/PostProcessing";
import SaveSelectedRunsToViewForm from "@/components/block/SaveSelectedRunsToViewForm";
//import RunActions from "@/components/block/RunActions";
import ListView from "../overlay/ListView.vue";
import LoadingOverlay from "../overlay/LoadingOverlay.vue";
import Badge from "../blocks/Badge.vue";
import { eventBus } from "@/event-bus";
import { producesPlottables } from "@/fileData";
import ComparePlots from "../blocks/ComparePlots.vue";


export default {
  name: 'ExperimentView',
  store: store,
  components: {
    ListView, LoadingOverlay, Badge,
    //RunActions,
    SaveSelectedRunsToViewForm, ComparePlots,
    //ButtonOverlay,
    //TableOverlayInfo,
    Multipane, MultipaneResizer,
    //PostProcessing
  },
  //store: store,
  //props: ['exp_name', 'status'],
  data() {
    return {
      selected: [],
      filterText: "",
      refreshDataInterval: null,
      view: null,

      //runIds: null,
      //selectedRunIdsMap: {},
      //runListRefreshTimeout: null,

      //processingList: false,
      //processingDeleteSelected: false,
      //processingDelete: {},

      //page: 1,
      //pageSize: 15
    };
  },
  computed: {
    isCompareEnabled() {
            return [...this.selected].some(
                run => {
                    let type1 = run.inputs.find(input => input.name == "Calculation_Type");
                    let type2 = run.inputs.find(input =>
                        ["EPOLYSCAT_Application_Module", "Application_Utility", "Application_Workflow"].indexOf(input.name) != -1
                    );

                    return run.status == "COMPLETED" &&
                            producesPlottables.indexOf(`${type1.value}/${type2.value}`.toUpperCase()) >= 0;
                }
            );
    },
    numberOfRunsSelected() {
                return this.selected.length;
    },
    isTutorials() {
                let splitPath = this.$route.path.split("/") || [""];

                return splitPath[splitPath.length - 1] == "tutorials";
    },
    viewId() {
                return (this.isTutorials) ? -1 :
                    ("viewId" in this.$route.params) ? this.$route.params.viewId : null;
    },

    //experimentId() {
    //  return this.$route.query.experimentId
    //},
    //viewId() {
    //  return this.$route.query.viewId
    //}, 
    //experiment() {
    //  return this.$store.getters["experiment/getExperiment"]({experimentId: this.experimentId});
    //},
    //view() {
    //  return this.$store.getters["view/getView"]({viewId: this.viewId});
    //},
    runs() {
         let runs = (this.view == null) ? this.$store.getters["run/getRuns"]() : this.view.runs;

         return runs.filter(run => run.name.includes(this.filterText)).map(run => {
                    run.isSelectable = !this.isTutorials;
                    return run;
         });
    },
//      if (this.processingList) {
//        return null;
//      } else {
//        return this.$store.getters["run/getRuns"]({
//          experimentId: this.experimentId,
//          viewId: this.viewId,
//          page: this.page,
//          pageSize: this.pageSize
//        });
//      }
//    }, 
    newRunLink() {  
        return (this.viewId != null) ? `/runs/new?viewId=${this.viewId}` : `/runs/new`;
    },
    loadingOverlayName() {
                return (this.viewId == null) ? "runs" : `viewRuns${this.viewId}`
    },
    pageName() {
                return (this.viewId == null) ? "All Runs" : (this.view == null) ? `Loading View...` : this.view.name
    },
    pageIsLoading() {
                return this.$store.getters['loading/getMessages'](this.loadingOverlayName).length > 0
    },
          
/*
    runsPagination() {
      return this.$store.getters["run/getRunsPagination"]({
        experimentId: this.experimentId,
        viewId: this.viewId,
        page: this.page,
        pageSize: this.pageSize
      });
    },
    selectedCount() {
      return this.selectedRunIds.length;
    },
    selectedRunIds() {
      const _selectedRunIds = [];
      for (let runId in this.selectedRunIdsMap) {
        if (this.selectedRunIdsMap[runId]) {
          _selectedRunIds.push(runId);
        }
      }

      return _selectedRunIds;
    }
*/
  },
  methods: {
    runLink({runId}) {
      let _link = "/runs/";
      if (runId) {
        _link += `${runId}?`;
      }

      if (this.experimentId) {
        _link += `experimentId=${this.experimentId}&`;
      }

      if (this.viewId) {
        _link += `viewId=${this.viewId}&`;
      }

      return _link;
    },
    async deleteAllSelectedRuns() {
      this.processingDeleteSelected = true;
      for (let i = 0; i < this.selectedRunIds.length; i++) {
        this.processingDelete = {...this.processingDelete, [this.selectedRunIds[i]]: true};
      }

      if (this.viewId) {
        try {
          await ViewService.removeRuns({viewId: this.viewId, runIds: this.selectedRunIds});
          await this.refreshList();
        } catch (e) {
          //TODO 
        }
      } else {
        try {
          await Promise.all(this.selectedRunIds.map(runId => {
            RunService.deleteRun({runId});
          }));
        } catch (e) {
          //TODO 
        }
      }

      this.selectedRunIdsMap = {};
      this.processingDeleteSelected = false;
      for (let i = 0; i < this.selectedRunIds.length; i++) {
        this.processingDelete = {...this.processingDelete, [this.selectedRunIds[i]]: false};
      }
    },
    async deleteRun({runId}) {
      this.processingDelete = {...this.processingDelete, [runId]: true};
      try {
        if (this.viewId) {
          await ViewService.removeRuns({viewId: this.viewId, runIds: [runId]});
        } else {
          await RunService.deleteRun({runId});
        }

        await this.refreshList();
      } catch (e) {
        //TODO 
      }
      this.processingDelete = {...this.processingDelete, [runId]: false};
    },
/*
    async refreshList({showOverlay = false} = {}) {
      if (showOverlay) this.processingList = true;

      await this.$store.dispatch("run/fetchRuns", {
        experimentId: this.experimentId,
        viewId: this.viewId,
        page: this.page,
        pageSize: this.pageSize
      });

      if (showOverlay) this.processingList = false;
    },
    sortStatus(status1, status2) {
                // orders both job_status and experiment status 
                const statuses = [
                    "UNSUBMITTED",
                    "CREATED",
                    "SUBMITTED",
                    "VALIDATED",
                    "QUEUED",
                    "SCHEDULED",
                    "LAUNCHED",
                    "EXECUTING",
                    "ACTIVE",
                    "SUSPENDED",
                    "CANCELING",
                    "CANCELED",
                    "COMPLETED",
                    "COMPLETE",
                    "NON_CRITICAL_FAIL",
                    "FAILED",
                    "UNKNOWN",
                ];

                return statuses.indexOf(status1.toUpperCase()) - statuses.indexOf(status2.toUpperCase());
    },
*/
    async refreshData(showLoading) {
         if (showLoading)        
             this.$store.commit("loading/START", { key: this.loadingOverlayName, message: "Fetching Runs" });
                                        
         if (this.viewId != null) {
                    try {
                        this.view = await this.$store.dispatch("view/fetchView", { viewId: this.viewId });
                    } catch (error) {
                        eventBus.$emit("error", { name: `Error while trying to fetch view with id: ${this.viewId}`, error });
                    }   
         } else {    
                    try {
                        await this.$store.dispatch("run/fetchRuns", {});
                    } catch (error) { 
                        eventBus.$emit("error", { name: `Error while trying to fetch the runs`, error }); 
                    }       
         }
            
         this.$store.commit("loading/STOP", { key: this.loadingOverlayName, message: "Fetching Runs" });
    },
    async removeFromView(runs) {
        if (this.view != null) {
            const removedRunIds = runs.map(run => run.id);

            this.$store.commit("loading/START", { key: this.loadingOverlayName, message: `Removing runs` });

            try {
                  await this.$store.dispatch("view/updateView", {
                    overide: true,
                    runIds: this.runs.map(run => run.id).filter(id => removedRunIds.indexOf(id) < 0),
                    viewId: this.viewId,
                    name: this.view.name
                  });
            } catch (error) {
                        eventBus.$emit("error", {
                            name: `Error while trying to remove runs from view with id: ${this.viewId}`,
                            error
                        });
            }

            this.$store.commit("loading/STOP", { key: this.loadingOverlayName, message: `Removing runs` });

            this.refreshData(true);
         }
    },
    async cloneRun(run) {
                await this.$router.push(`/runs/new?clonedFrom=${run.id}`);
                this.$router.go()
    },
    async deleteRuns(runs) {
        const runOrRuns = runs.length == 1 ? "run" : "runs";

        const deleteAll = await this.$bvModal.msgBoxConfirm(
               `Are you sure? Deleting ${runs.length} ${runOrRuns}: ${
                   runs
                       .map(run => `"${run.name}"`)
                       .slice(0, 4)
                       .join(",")
                } ${runs.length > 4 ? "..." : "" }`, {
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
                   `Do want to delete all files associated with the ${runOrRuns}?`, {
                       title: 'Please Confirm',
                       okVariant: 'danger',
                       okTitle: 'Yes',
                       cancelTitle: 'No',
                       footerClass: 'p-2',
                       hideHeaderClose: false,
                       centered: true
                    }
               )

               for (const run of runs) {
                        this.$store.commit("loading/START", { key: this.loadingOverlayName, message: `Deleting "${run.name}"` });

                  try {
                            await this.$store.dispatch("run/deleteRun", {
                                runId: run.id,
                                deleteAssociated
                            });

                   } catch (error) {
                        eventBus.$emit("error", {
                           name: `Error while trying to delete run with id: ${run.id}`,
                                error
                        });
                   }

                   this.$store.commit("loading/STOP", { key: this.loadingOverlayName, message: `Deleting "${run.name}"` });
              }

              this.refreshData(true);
        }
    },
    async submitRun(run) {
        this.$store.commit("loading/START", { key: this.loadingOverlayName, message: `Submitting "${run.name}"` });

           try {
                 run = await this.$store.dispatch("run/submitRun", {
                        runId: run.id
                 });
           } catch (error) {
                    eventBus.$emit("error", { name: `Error while trying to submit the run`, error });
           }

           this.$store.commit("loading/STOP", { key: this.loadingOverlayName, message: `Submitting "${run.name}"` });
    },
    async resubmitRun(run) {
            this.$store.commit("loading/START", { key: this.loadingOverlayName, message: `Resubmitting "${run.name}"` });

                try {
                    run = await this.$store.dispatch("run/submitRun", {
                        runId: run.id
                    });
                } catch (error) {
                    eventBus.$emit("error", { name: `Error while trying to submit the run`, error });
                }

                this.$store.commit("loading/STOP", { key: this.loadingOverlayName, message: `Resubmitting "${run.name}"` });
    },
            sortStatus(status1, status2) {
                // orders both job_status and experiment status
                const statuses = [
                    "UNSUBMITTED",
                    "CREATED",
                    "SUBMITTED",
                    "VALIDATED",
                    "QUEUED",
                    "SCHEDULED",
                    "LAUNCHED",
                    "EXECUTING",
                    "ACTIVE",
                    "SUSPENDED",
                    "CANCELING",
                    "CANCELED",
                    "COMPLETED",
                    "COMPLETE",
                    "NON_CRITICAL_FAIL",
                    "FAILED",
                    "UNKNOWN",
                ];

                return statuses.indexOf(status1.toUpperCase()) - statuses.indexOf(status2.toUpperCase());
            }
  },
  mounted() {
      if (this.viewId != null)
                this.refreshData(true);

      this.refreshDataInterval = setInterval(() => {
                this.refreshData(false);
            }, 45000);
      },
      beforeDestroy() {
      clearInterval(this.refreshDataInterval);
    }

}
</script>

<style scoped>
    .filter {
        margin-bottom: 20px;
        margin-top: -5px;
        width: 93%;
    }

    .button_row > * {
        margin: 0 10px;

    }

    .multipane .multipane-resizer {
      margin: 0;
      left: 0;
      position: relative;
    }

    .multipane-resizer {
        width: 0;
        border-left: 1px solid #ddd;
        pointer-events: none;
    }

</style>
