<template>
  <div class="w-100 h-100 bg-light p-2">
    <div class="w-100 h-100 bg-white overflow-auto p-3">
      <b-breadcrumb>
        <b-breadcrumb-item to="/views">Views</b-breadcrumb-item>
      </b-breadcrumb>
      <div class="w-100 mb-3">
        <div class="d-inline" style="font-weight: 400; font-size: 19px;">Views</div>
      </div>
      <div class="w-100 text-right mb-2">
        <!--        <router-link to="/create-view" v-slot="{ href, route, navigate, isActive,isExactActive }">-->
        <!--          <b-button variant="primary" size="sm" tag="a" :class="{active: isExactActive}" :href="href"-->
        <!--                    @click="navigate">-->
        <!--            Create new view-->
        <!--          </b-button>-->
        <!--        </router-link>-->
        <button-overlay :show="processingDeleteSelected">
          <b-button variant="outline-primary" size="sm" v-if="selectedCount > 0" class="ml-2"
                    v-b-modal="'modal-confirmation-delete-selected-views'">
            <b-icon icon="trash"></b-icon>
            Delete selected ({{ selectedCount }}) view{{ selectedCount > 1 ? 's' : '' }}
          </b-button>
        </button-overlay>
        <b-modal id="modal-confirmation-delete-selected-views" title="Delete Confirmation" ok-title="Delete"
                 v-on:ok="deleteAllSelectedViews">
          <p class="my-4">Are you sure that you want to delete {{ selectedCount }} selected views? </p>
        </b-modal>
      </div>
      <div class="w-100 overflow-auto">
        <table-overlay-info :data="views" :rows="5" :columns="4" empty-label="No views available.">
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
              <b-tr v-for="view in views" :key="view.viewId">
                <b-td>
                  <b-form-checkbox v-model="selectedViewIdsMap[view.viewId]"/>
                </b-td>
                <b-td>
                  <router-link :to="`/runs/?viewId=${view.viewId}`"
                               v-slot="{ href, route, navigate, isActive,isExactActive }">
                    <b-link :class="{active: isExactActive}" :href="href" @click="navigate">
                      <div class="overflow-auto" style="max-width: 200px;">{{ view.name }}</div>
                    </b-link>
                  </router-link>
                </b-td>
                <b-td>
                  <div class="overflow-auto" style="max-width: 200px;">{{ view.description }}</div>
                </b-td>
                <b-td>{{ view.updated }}</b-td>
                <b-td class="text-center" style="min-width: 110px; max-width: 110px;">
                  <router-link :to="`/runs/?viewId=${view.viewId}`"
                               v-slot="{ href, route, navigate, isActive,isExactActive }">
                    <b-link :href="href" @click="navigate"> {{ view.activeRunCount }}</b-link>
                  </router-link>
                </b-td>
                <b-td style="min-width: 100px; max-width: 100px;">
                  <button-overlay :show="processingDelete[view.viewId]">
                    <b-button variant="link" size="sm" class="ml-2" v-b-tooltip.hover.auto title="Delete"
                              v-b-modal="`modal-confirmation-delete-${view.viewId}`">
                      <b-icon icon="trash"></b-icon>
                    </b-button>
                  </button-overlay>
                  <b-modal :id="`modal-confirmation-delete-${view.viewId}`" title="Delete Confirmation"
                           ok-title="Delete"
                           v-on:ok="deleteView(view)">
                    <p class="my-4">Are you sure that you want to delete the view "{{ view.name }}"? </p>
                  </b-modal>
                </b-td>
              </b-tr>
            </b-tbody>
          </b-table-simple>
        </table-overlay-info>
      </div>
      <b-pagination
          v-model="page"
          :total-rows="viewsPagination.total"
          :per-page="pageSize"
      ></b-pagination>
    </div>
  </div>
</template>

<script>
import store from "@/store";
import TableOverlayInfo from "@/components/overlay/table-overlay-info";
import ButtonOverlay from "@/components/overlay/button-overlay";
import {ViewService} from "@/service/epolyscat-service";

export default {
  name: "Views",
  components: {TableOverlayInfo, ButtonOverlay},
  store: store,
  data() {
    return {
      page: 1,
      pageSize: 10,

      selectedViewIdsMap: {},
      processingDelete: {},
      processingDeleteSelected: false
    }
  },
  computed: {
    views() {
      return this.$store.getters["view/getViews"]({page: this.page, pageSize: this.pageSize});
    },
    viewsPagination() {
      return this.$store.getters["view/getViewsPagination"]({page: this.page, pageSize: this.pageSize});
    },
    selectedCount() {
      return this.selectedViewIds.length;
    },
    selectedViewIds() {
      const _selectedViewIds = [];
      for (let viewId in this.selectedViewIdsMap) {
        if (this.selectedViewIdsMap[viewId]) {
          _selectedViewIds.push(viewId);
        }
      }

      return _selectedViewIds;
    }
  },
  methods: {
    refreshData() {
      this.$store.dispatch("view/fetchViews", {page: this.page, pageSize: this.pageSize});
    },
    deleteAllSelectedViews() {
      this.processingDeleteSelected = true;
      for (let i = 0; i < this.selectedViewIds.length; i++) {
        this.deleteView({viewId: this.selectedViewIds[i]});
      }
      this.selectedViewIdsMap = {};
      this.processingDeleteSelected = false;
    },
    async deleteView({viewId}) {
      this.processingDelete = {...this.processingDelete, [viewId]: true};
      try {
        await ViewService.deleteView({viewId});
        this.refreshData();
      } catch (e) {
        //TODO
      }
      this.processingDelete = {...this.processingDelete, [viewId]: false};
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
