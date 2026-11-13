<template>
  <div class="portal-page">
    <main class="portal-page-content">
      <header class="portal-page-header"><h1 class="portal-page-title">Documentation</h1></header>
      <p>Read the ePolyScat manual for input records, commands, and scientific examples.</p>
      <p><a href="https://rrlucchese.github.io/">Open the ePolyScat manual</a></p>
      <iframe class="portal-documentation-frame" title="ePolyScat manual" src="https://rrlucchese.github.io/" />
    </main>
  </div>
</template>

<script>
export default {
  data() {
        return {
          embeddedUrl: 'https://www.example.com' // Initial URL
        };
      }
};
  /*
  },
  methods: {
    async fetchHtml() {
      try {
        const response = await fetch('https://rrlucchese.github.io/');
        this.externalHtml = await response.text();
      } catch (error) {
        console.error('Error fetching HTML:', error);
      }
    }
  }
  */
/*
export default {
  name: 'Home',
  store: store,
  //components: { ExperimentCard },
  data() {
    return {
      numberOfExperiments: 5,
      experimentStatistics: {
        experimentCount: null,
        runCount: null
      }
    }
  },
  computed: {
    experiments() {
      return this.$store.getters["experiment/getExperiments"]({ pageSize: this.numberOfExperiments });
    }
  },
  methods: {
    async refreshData() {
      this.$store.dispatch("experiment/fetchExperiments", { pageSize: this.numberOfExperiments });
      this.$store.dispatch("run/fetchRuns");

      this.experimentStatistics = await ExperimentService.fetchExperimentStatistics();
    }
  },
  mounted() {
    this.refreshData();
  }
}

*/
</script>

<style scoped>
.stats-card {
  background: #2265970D;
}

.exp-card {
  width: 300px !important;
  height: 180px !important;
}

.create-exp-card {
  border: 2px dashed #1C6296;
  box-sizing: border-box;
  border-radius: 5px;
  font-weight: bold;
  line-height: 70px;
}

.grid-container {
  column-gap: 20px !important;
  row-gap: 40px !important;
  display: grid;
  grid-auto-columns: auto;
  grid-template-columns: repeat( auto-fit, 300px);
}
</style>
