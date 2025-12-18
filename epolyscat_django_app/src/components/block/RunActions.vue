<template>
  <div class="d-inline">
    <router-link :to="runCloneLink(run)"
                 v-slot="{ href, route, navigate, isActive,isExactActive }">
      <b-button :variant="variant" size="sm" @click="navigate" v-b-tooltip.hover.auto
                :title="showButtonText? '' : 'Clone'" class="mb-2">
        <b-icon icon="clipboard"></b-icon>
        <span class="pl-2" v-if="showButtonText">Clone</span>
      </b-button>
    </router-link>

    <template v-if="!view || view.type !== 'tutorial'">
      <router-link :to="runResubmitLink(run)"
                   v-slot="{ href, route, navigate, isActive,isExactActive }">

        <template v-if="run">
          <b-button :variant="variant" size="sm" @click="navigate" v-b-tooltip.hover.auto
                    :title="showButtonText? '' : 'Submit'" class="mb-2"
                    v-if="run.canSubmit">
            <b-icon icon="arrow-bar-up"></b-icon>
            <span class="pl-2" v-if="showButtonText">Submit to server</span>
          </b-button>

          <b-button :variant="variant" size="sm" @click="navigate" v-b-tooltip.hover.auto
                    :title="showButtonText? '' : 'Resubmit'" class="mb-2"
                    v-else :disabled="!run.canResubmit">
            <b-icon icon="arrow90deg-right"></b-icon>
            <span class="pl-2" v-if="showButtonText">Resubmit to server</span>
          </b-button>
        </template>

      </router-link>
    </template>

    <template v-if="!view || !view.readonly">
      <button-overlay :show="processingDelete[run.runId]">
        <b-button :variant="variant" size="sm" class="ml-1 mb-2" v-b-tooltip.hover.auto
                  :title="showButtonText ? '' : view ? 'Remove from selection' : 'Delete'"
                  v-on:click="onDeleteClick">

          <b-icon v-if="view" icon="x"></b-icon>
          <b-icon v-else icon="trash"></b-icon>

          <span class="pl-2" v-if="showButtonText">
            <template v-if="view">Remove</template>
            <template v-else>Destroy</template>
          </span>
        </b-button>
      </button-overlay>
      <b-modal :id="`modal-confirmation-delete-${run.runId}`" title="Delete Confirmation"
               ok-title="Delete"
               v-on:ok="deleteRun(run)">
        <p class="my-4">Are you sure that you want to delete the run "{{ run.name }}"? </p>
      </b-modal>
    </template>
  </div>
</template>

<script>
import {RunService, ViewService} from "@/service/epolyscat-service";
import ButtonOverlay from "@/components/overlay/button-overlay";

export default {
  name: "RunActions",
  components: {ButtonOverlay},
  props: {
    "run": {}, "view": {}, "experiment": {},
    "variant": {default: "link"},
    "showButtonText": {default: false}
  },
  data() {
    return {
      processingDelete: false,
    }
  },
  methods: {
    runCloneLink({runId}) {
      let _link = "/create-run?";
      if (runId) {
        _link += `cloneRunId=${runId}&`;
      }

      if (this.experiment) {
        _link += `experimentId=${this.experiment.experimentId}&`;
      }

      if (this.view) {
        _link += `viewId=${this.view.viewId}&`;
      }

      return _link;
    },
    runResubmitLink({runId}) {
      let _link = "/resubmit-run?";
      if (runId) {
        _link += `runId=${runId}&`;
      }

      if (this.experiment) {
        _link += `experimentId=${this.experiment.experimentId}&`;
      }

      if (this.view) {
        _link += `viewId=${this.view.viewId}&`;
      }

      return _link;
    },
    onDeleteClick() {
      if (this.view) {
        this.deleteRun(this.run);
      } else {
        this.$bvModal.show(`modal-confirmation-delete-${this.run.runId}`);
      }
    },
    async deleteRun({runId}) {
      this.processingDelete = {...this.processingDelete, [runId]: true};
      try {
        if (this.view) {
          await ViewService.removeRuns({viewId: this.view.viewId, runIds: [runId]});
        } else {
          await RunService.deleteRun({runId});
        }

        await this.$emit("delete", runId);
      } catch (e) {
        //TODO
      }
      this.processingDelete = {...this.processingDelete, [runId]: false};
    },
  }
}
</script>

<style scoped>

</style>
