<template>
  <div>
    <div>
      <label>View name*</label>
      <b-form-input
          id="root"
          size="sm"
          v-model="selectedViewName"
          :state="inputState.selectedViewName"
          list="root-list"
          autocomplete="off"
      />
      <b-form-datalist id="root-list" :options="viewNames"/>
      <b-form-invalid-feedback>
        The root of the run cannot be empty
      </b-form-invalid-feedback>
    </div>
    <div class="w-100 mt-2">
      <b-button class="w-100" variant="primary" size="sm" v-on:click="onSaveClick">Save</b-button>
    </div>
  </div>
</template>

<script>
import store from "@/store";

export default {
  name: "SaveSelectedRunsToViewForm",
  store: store,
  props: ["runIds"],
  data() {
    return {
      selectedViewName: null,

      inputFieldsList: ["selectedViewName"]
    }
  },
  computed: {
    inputState() {
      return {
        selectedViewName: this.selectedViewName === null ? null : this.isValid.selectedViewName
      }
    },
    isValid() {
      return {
        selectedViewName: !!this.selectedViewName && this.selectedViewName.length > 0
      }
    },
    isFormValid() {
      let _isFormValid = true;
      for (let i = 0; i < this.inputFieldsList.length; i++) {
        _isFormValid = _isFormValid && this.isValid[this.inputFieldsList[i]];
      }

      return _isFormValid;
    },
    views() {
      return this.$store.getters["view/getViews"]();
    },
    writableViews() {
      if (this.views) {
        return this.views.filter(({readonly}) => !readonly);
      }

      return [];
    },
    viewNames() {
      return this.writableViews.map(({name}) => name);
    }
  },
  methods: {
    makeFormVisited() {
      for (let i = 0; i < this.inputFieldsList.length; i++) {
        if (this[this.inputFieldsList[i]] === null) {
          this[this.inputFieldsList[i]] = "";
        }
      }
    },
    async onSaveClick() {
      this.makeFormVisited();
      if (this.isFormValid) {
        const viewIndex = this.viewNames.indexOf(this.selectedViewName);
        let view;
        if (viewIndex < 0) {
          view = await this.$store.dispatch("view/createView", {
            name: this.selectedViewName,
            runIds: this.runIds
          });
        } else {
          view = await this.$store.dispatch("view/updateView", {
            viewId: this.writableViews[viewIndex].viewId,
            name: this.selectedViewName,
            runIds: this.runIds
          });
        }

        if (view) {
          await this.$store.dispatch("view/fetchViews");
          this.$emit('close');
        }
      }
    }
  },
  mounted() {
    this.$store.dispatch("view/fetchViews");
  }
}
</script>

<style scoped>

</style>