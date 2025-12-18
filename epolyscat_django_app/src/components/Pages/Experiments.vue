<template>
  <div class="w-100 h-100 bg-light p-2">
    <div class="w-100 h-100 bg-white overflow-auto p-3">
      <b-breadcrumb>
        <b-breadcrumb-item to="/experiments">Experiments</b-breadcrumb-item>
      </b-breadcrumb>
      <div class="w-100 mb-3">
        <div class="d-inline" style="font-weight: 400; font-size: 19px;">Experiments</div>
      </div>
      <div class="w-100 text-right mb-2">
        <router-link to="/create-experiment" v-slot="{ href, route, navigate, isActive,isExactActive }">
          <b-button variant="primary" size="sm" tag="a" :class="{active: isExactActive}" :href="href"
                    @click="navigate">
            Create new experiment
          </b-button>
        </router-link>
        <button-overlay :show="processingDeleteSelected">
          <b-button variant="outline-primary" size="sm" v-if="selectedCount > 0" class="ml-2"
                    v-b-modal="'modal-confirmation-delete-selected-experiments'">
            <b-icon icon="trash"></b-icon>
            Delete selected ({{ selectedCount }}) experiment{{ selectedCount > 1 ? 's' : '' }}
          </b-button>
        </button-overlay>
        <b-modal id="modal-confirmation-delete-selected-experiments" title="Delete Confirmation" ok-title="Delete"
                 v-on:ok="deleteAllSelectedExperiments">
          <p class="my-4">Are you sure that you want to delete {{ selectedCount }} selected experiments? </p>
        </b-modal>
      </div>
      <div class="w-100 overflow-auto">
        <table-overlay-info :data="experiments" :rows="5" :columns="4" empty-label="No experiments available.">
          <b-table-simple>
            <b-thead>
              <b-tr>
                <b-th></b-th>
                <b-th>Name</b-th>
                <b-th>Description</b-th>
                <b-th>Last edited</b-th>
                <b-th>Active runs</b-th>
                <b-th>Actions</b-th>
              </b-tr>
            </b-thead>
            <b-tbody>
              <b-tr v-for="experiment in experiments" :key="experiment.experimentId">
                <b-td>
                  <b-form-checkbox v-model="selectedExperimentIdsMap[experiment.experimentId]"/>
                </b-td>
                <b-td>
                  <router-link :to="`/runs/?experimentId=${experiment.experimentId}`"
                               v-slot="{ href, route, navigate, isActive,isExactActive }">
                    <b-link :class="{active: isExactActive}" :href="href" @click="navigate">
                      <div class="overflow-auto" style="max-width: 200px;">{{ experiment.name }}</div>
                    </b-link>
                  </router-link>
                </b-td>
                <b-td>
                  <div class="overflow-auto" style="max-width: 200px;">{{ experiment.description }}</div>
                </b-td>
                <b-td>{{ experiment.updated }}</b-td>
                <b-td class="text-center" style="min-width: 110px; max-width: 110px;">
                  <router-link :to="`/runs/?experimentId=${experiment.experimentId}`"
                               v-slot="{ href, route, navigate, isActive,isExactActive }">
                    <b-link :href="href" @click="navigate"> {{ experiment.activeRunCount }}</b-link>
                  </router-link>
                </b-td>
                <b-td style="min-width: 100px; max-width: 100px;">
                  <button-overlay :show="processingDelete[experiment.experimentId]">
                    <b-button variant="link" size="sm" class="ml-2" v-b-tooltip.hover.auto title="Delete"
                              v-b-modal="`modal-confirmation-delete-${experiment.experimentId}`">
                      <b-icon icon="trash"></b-icon>
                    </b-button>
                  </button-overlay>
                  <b-modal :id="`modal-confirmation-delete-${experiment.experimentId}`" title="Delete Confirmation"
                           ok-title="Delete"
                           v-on:ok="deleteExperiment(experiment)">
                    <p class="my-4">Are you sure that you want to delete the experiment "{{ experiment.name }}"? </p>
                  </b-modal>
                </b-td>
              </b-tr>
            </b-tbody>
          </b-table-simple>
        </table-overlay-info>
      </div>
      <b-pagination
          v-model="page"
          :total-rows="experimentsPagination.total"
          :per-page="pageSize"
      ></b-pagination>
    </div>
  </div>
</template>

<script>
import store from "@/store";
import TableOverlayInfo from "@/components/overlay/table-overlay-info";
import ButtonOverlay from "@/components/overlay/button-overlay";
import {ExperimentService} from "@/service/epolyscat-service";

export default {
  name: "Experiments",
  components: {TableOverlayInfo, ButtonOverlay},
  store: store,
  data() {
    return {
      page: 1,
      pageSize: 10,

      selectedExperimentIdsMap: {},
      processingDelete: {},
      processingDeleteSelected: false
    }
  },
  computed: {
    experiments() {
      return this.$store.getters["experiment/getExperiments"]({page: this.page, pageSize: this.pageSize});
    },
    experimentsPagination() {
      return this.$store.getters["experiment/getExperimentsPagination"]({page: this.page, pageSize: this.pageSize});
    },
    selectedCount() {
      return this.selectedExperimentIds.length;
    },
    selectedExperimentIds() {
      const _selectedExperimentIds = [];
      for (let experimentId in this.selectedExperimentIdsMap) {
        if (this.selectedExperimentIdsMap[experimentId]) {
          _selectedExperimentIds.push(experimentId);
        }
      }

      return _selectedExperimentIds;
    }
  },
  methods: {
    refreshData() {
      this.$store.dispatch("experiment/fetchExperiments", {page: this.page, pageSize: this.pageSize});
    },
    deleteAllSelectedExperiments() {
      this.processingDeleteSelected = true;
      for (let i = 0; i < this.selectedExperimentIds.length; i++) {
        this.deleteExperiment({experimentId: this.selectedExperimentIds[i]});
      }
      this.selectedExperimentIdsMap = {};
      this.processingDeleteSelected = false;
    },
    async deleteExperiment({experimentId}) {
      this.processingDelete = {...this.processingDelete, [experimentId]: true};
      try {
        await ExperimentService.deleteExperiment({experimentId});
        this.refreshData();
      } catch (e) {
        //TODO
      }
      this.processingDelete = {...this.processingDelete, [experimentId]: false};
    },
  },
  watch: {
    page() {
      this.refreshData();
    },
    pageSize() {
      this.refreshData();
    }
  },
  mounted() {
    this.refreshData();
  }
}
</script>

<style scoped>

</style>
