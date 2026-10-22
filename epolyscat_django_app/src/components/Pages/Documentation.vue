<template>
  <div class="w-100 h-100 overflow-auto">
    <div class="w-100 d-flex p-5">
      <div class="flex-fill">
        <div class="text-primary" style="font-size: 56px;">ePolyScat Documentation</div>
        <div>
          ePolyScat is an electron–molecule scattering and photoionization
          suite for computing photoelectron spectra, differential cross sections,
          resonances, and angular distributions. This portal lets you create,
          submit, and analyze ePS runs on modern HPC resources.
        </div>
      </div>
    </div>
    <div class="w-100">
      <iframe
        src="https://epolyscat.droppages.com/IntroREADME2/"
        title="ePolyScat documentation"
        width="100%"
        height="700"
        frameborder="0"
      ></iframe>
    </div>
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
        const response = await fetch('https://epolyscat.droppages.com/front');
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
