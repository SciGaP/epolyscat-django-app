<template>
  <div class="w-100 h-100 bg-light p-2">
    <multipane class="w-100 h-100" layout="vertical">
      <div class="h-100 bg-white overflow-auto" style="flex-grow: 1">
        <div class="w-100 p-3">
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
              <router-link :to="`/create-run?experimentId=${experimentId}`"
                           v-if="!viewId && runIds && runIds.length === 0"
                           v-slot="{ href, route, navigate, isActive,isExactActive }">
                <b-button variant="primary" size="sm" tag="a" :class="{active: isExactActive}" :href="href"
                          @click="navigate">
                  Add a new run
                </b-button>
              </router-link>

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
        </div>
        <div class="w-100 p-3 overflow-auto">
          <table-overlay-info :data="runs" :rows="5" :columns="5" empty-label="No runs available.">
            <b-table-simple>
              <b-thead>
                <b-tr>
                  <b-th></b-th>
                  <b-th>Name</b-th>
                  <b-th v-if="!view || view.type !== 'tutorial'">Status</b-th>
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
            :per-page="pageSize"
        ></b-pagination>
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
      </div>
      <multipane-resizer v-if="selectedCount > 0"/>
      <div class="h-100 m-1" v-if="selectedCount > 0">
        <PostProcessing :selectedRunIds="selectedRunIds"/>
      </div>
    </multipane>
  </div>
</template>

<script>
import store from "@/store";
import TableOverlayInfo from "@/components/overlay/table-overlay-info";
import {RunService, ViewService} from "@/service/trecx-service";
import ButtonOverlay from "@/components/overlay/button-overlay";
import {Multipane, MultipaneResizer} from 'vue-multipane';
import PostProcessing from "@/components/block/PostProcessing";
import SaveSelectedRunsToViewForm from "@/components/block/SaveSelectedRunsToViewForm";
import RunActions from "@/components/block/RunActions";

export default {
  name: 'ExperimentView',
  components: {
    RunActions,
    SaveSelectedRunsToViewForm,
    ButtonOverlay,
    TableOverlayInfo,
    Multipane,
    MultipaneResizer,
    PostProcessing
  },
  store: store,
  props: ['exp_name', 'status'],
  data() {
    return {
      runIds: null,
      selectedRunIdsMap: {},
      runListRefreshTimeout: null,

      processingList: false,
      processingDeleteSelected: false,
      processingDelete: {},

      page: 1,
      pageSize: 15
    }
  },
  computed: {
    experimentId() {
      return this.$route.query.experimentId
    },
    viewId() {
      return this.$route.query.viewId
    },
    experiment() {
      return this.$store.getters["experiment/getExperiment"]({experimentId: this.experimentId});
    },
    view() {
      return this.$store.getters["view/getView"]({viewId: this.viewId});
    },
    runs() {
      if (this.processingList) {
        return null;
      } else {
        return this.$store.getters["run/getRuns"]({
          experimentId: this.experimentId,
          viewId: this.viewId,
          page: this.page,
          pageSize: this.pageSize
        });
      }
    },
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
    async refreshData() {
      if (this.experimentId) {
        await this.$store.dispatch("experiment/fetchExperiment", {
          experimentId: this.experimentId
        });
      }

      if (this.viewId) {
        await this.$store.dispatch("view/fetchView", {
          viewId: this.viewId
        });
      }

      await this.refreshList({showOverlay: true});
    }
  },
  watch: {
    page() {
      this.refreshList({showOverlay: true});
    },
    pageSize() {
      this.refreshList({showOverlay: true});
    },
    experimentId() {
      this.selectedRunIdsMap = {};
      this.refreshList({showOverlay: true});
    },
    viewId() {
      this.selectedRunIdsMap = {};
      this.refreshList({showOverlay: true});
    },
    runs() {
      if (this.runs) {
        this.runIds = this.runs.sort(function (a, b) {
          return new Date(b.created) - new Date(a.created);
        }).map(({runId}) => runId);
      }
    },
    selectedRunIds() {
      this.$emit('on-select-runs', this.selectedRunIds);
    }
  },
  mounted() {
    this.$emit('on-select-runs', this.selectedRuns);
    this.refreshData();

    this.runListRefreshTimeout = setInterval(() => {
      this.refreshList();
    }, 15000);
  },
  beforeDestroy() {
    clearInterval(this.runListRefreshTimeout);
  }
}
</script>

<style scoped>
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
</style>
