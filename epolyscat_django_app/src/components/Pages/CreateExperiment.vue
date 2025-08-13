<template>
  <div class="w-100 h-100 bg-light p-2">
    <div class="w-100 h-100 bg-white d-flex flex-column overflow-auto p-3">
      <b-breadcrumb>
        <b-breadcrumb-item to="/experiments">Experiments</b-breadcrumb-item>
        <b-breadcrumb-item to="/create-experiment">New</b-breadcrumb-item>
      </b-breadcrumb>

      <div class="w-100">
        <div class="d-inline" style="font-weight: 400; font-size: 19px;">Create a new experiment</div>
      </div>

      <div class="w-100 d-flex flex-row">
        <div class="p-2 d-flex flex-row flex-fill">
          <div style="min-width: 200px;"><label for="name">Name of the experiment *</label></div>
          <div>
            <b-form-input id="name" size="sm" v-model="name" :state="inputState.name"/>
            <b-form-invalid-feedback>
              The name of the experiment cannot be empty
            </b-form-invalid-feedback>
          </div>
        </div>
        <div class="p-2 d-flex flex-row flex-fill">
          <div style="min-width: 200px;"><label for="description">Description</label></div>
          <div>
            <b-form-input id="name" size="sm" v-model="description" :state="inputState.description"/>
          </div>
        </div>
      </div>

      <div class="w-100 p-2">
        <b-button variant="primary" class="w-100" v-on:click="onSubmit">Submit</b-button>
      </div>

    </div>
  </div>
</template>

<script>
import store from "@/store";

export default {
  name: "CreateExperiment",
  store: store,
  data() {
    return {
      name: null,
      description: null,

      inputFieldsList: ["name", "description"]
    }
  },
  computed: {
    inputState() {
      return {
        name: this.name === null ? null : this.isValid.name,
        description: this.description === null ? null : this.isValid.description
      }
    },
    isValid() {
      return {
        name: !!this.name && this.name.length >= 1,
        description: true
      }
    },
    isFormValid() {
      let _isFormValid = true;
      for (let i = 0; i < this.inputFieldsList.length; i++) {
        _isFormValid = _isFormValid && this.isValid[this.inputFieldsList[i]];
      }

      return _isFormValid;
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
    async onSubmit() {
      this.makeFormVisited();
      if (this.isFormValid) {
        const experiment = await this.$store.dispatch("experiment/createExperiment", {
          name: this.name,
          description: this.description
        });
        const {experimentId} = experiment;
        this.$router.history.push(`/runs/?experimentId=${experimentId}`);

        // For the left navigation to be updated
        this.$store.dispatch("experiment/fetchExperiments");
      }
    },
  }
}
</script>

<style scoped>

</style>