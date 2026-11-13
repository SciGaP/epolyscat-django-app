<template>
  <div class="portal-page">
    <div class="portal-page-content">
      <b-breadcrumb>
        <b-breadcrumb-item to="/experiments">Experiments</b-breadcrumb-item>
        <b-breadcrumb-item to="/create-experiment">New</b-breadcrumb-item>
      </b-breadcrumb>

      <header class="portal-page-header"><h1 class="portal-page-title">New Experiment</h1></header>

      <div class="portal-experiment-fields">
        <div class="portal-experiment-field">
          <div><label for="name">Name of the experiment *</label></div>
          <div>
            <b-form-input id="name" v-model="name" :state="inputState.name"/>
            <b-form-invalid-feedback>
              The name of the experiment cannot be empty
            </b-form-invalid-feedback>
          </div>
        </div>
        <div class="portal-experiment-field">
          <div><label for="description">Description</label></div>
          <div>
            <b-form-input id="description" v-model="description" :state="inputState.description"/>
          </div>
        </div>
      </div>

      <div class="portal-page-actions">
        <b-button variant="primary" v-on:click="onSubmit">Create Experiment</b-button>
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