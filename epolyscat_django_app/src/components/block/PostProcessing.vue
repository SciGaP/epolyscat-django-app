<template>
  <div class="w-100 h-100 p-2 bg-white overflow-auto" style="max-width: 400px;">
    <div class="w-100 p-2">
      <div style="font-weight: 400; font-size: 19px;">Post Processing</div>
    </div>

    <template v-if="selectedRunIds && selectedRunIds.length > 1">
      <div class="w-100 p-2">
        <b-button size="sm" v-for="outputOption in outputOptions" :key="outputOption" class="mb-2"
                  :variant="selectedOutputOptions[outputOption] ? 'primary': 'outline-secondary'"
                  v-on:click="selectedOutputOptions[outputOption] = !selectedOutputOptions[outputOption]">
          {{ outputOption }}
        </b-button>
      </div>
      <div class="w-100 p-2">
        <b-button size="sm" variant="outline-primary" v-on:click="selectAllOutputOptions" class="mb-2">
          Display All
        </b-button>
        <b-button size="sm" variant="outline-primary" v-on:click="unselectAllOutputOptions" class="mb-2">
          Unselect All
        </b-button>
      </div>
    </template>

    <Errors :errors="errors"/>

    <div class="w-100 p-2" v-if="hasPlotSelected">
      <b-overlay :show="processingPlot" rounded spinner-small>
        <img class="w-100" :src="plotImageUrl" v-if="plotImageUrl"/>
        <div class="text-center p-3 bg-light" v-else>
          (Plot)
        </div>
      </b-overlay>
    </div>

    <b-overlay :show="processingPlot || processingListInputs || processingInputDifference" rounded spinner-small>
      <div class="w-100 d-flex flex-row">
        <div class="p-2">
          <label for="expec">Data to plot</label>
          <b-overlay :show="!expecItemList" rounded spinner-small class="p-2">
            <b-form-select id="expec" size="sm" v-model="expectationValue" :options="expecItemList"/>
          </b-overlay>
        </div>
        <div class="p-2">
          <label for="params">With Parameters Below</label>
          <b-overlay :show="!parameterList" rounded spinner-small class="p-2">
            <b-form-select id="params" size="sm" v-model="selectedParameters" :options="parameterList"/>
          </b-overlay>
        </div>
      </div>
      <div class="w-100 d-flex flex-row">
        <div class="p-2 d-flex flex-row">
          <div style="min-width: 45px;"><label for="xAxis">X axis</label></div>
          <div>
            <b-form-input id="xAxis" size="sm" v-model="xAxis"/>
          </div>
        </div>
        <div class="p-2 d-flex flex-row">
          <div style="min-width: 45px;"><label for="yAxis">Y axis</label></div>
          <div>
            <b-form-input id="yAxis" size="sm" v-model="yAxis"/>
          </div>
        </div>
      </div>
      <div class="w-100 d-flex flex-row">
        <div class="p-2 d-flex flex-row" style="flex: 1;">
          <div style="min-width: 45px;"><label for="flags">Flags</label></div>
          <div>
            <b-form-input id="flags" size="sm" v-model="flags"/>
          </div>
        </div>
        <div class="p-2 d-flex flex-row" style="flex: 1;">
          <b-overlay :show="!expecItemList || !parameterList" rounded spinner-small>
            <b-button size="sm" variant="primary" class="w-100" v-on:click="plotRuns">Plot</b-button>
          </b-overlay>
        </div>
      </div>
    </b-overlay>

    <div class="w-100 p-2" v-if="hasListInputsSelected">
      <label>List Inputs</label>
      <b-overlay :show="processingListInputs" rounded spinner-small>
        <div class="w-100 p-2 overflow-auto d-inline-block console-output" v-if="listInputsOutput !== null">
          {{ listInputsOutput }}
        </div>
        <div class="text-center p-3 bg-light" v-else>
          (List Inputs)
        </div>
      </b-overlay>
    </div>

    <div class="w-100 p-2" v-if="hasInputDifferenceSelected">
      <label>Input Differences</label>
      <b-overlay :show="processingInputDifference" rounded spinner-small>
        <div class="w-100 p-2 overflow-auto d-inline-block console-output" v-if="inputDifferenceOutput !== null">
          {{ inputDifferenceOutput }}
        </div>
        <div class="text-center p-3 bg-light" v-else>
          (Input Differences)
        </div>
      </b-overlay>
    </div>
  </div>
</template>

<script>
import {PlotService} from "@/service/epolyscat-service";
import store from "@/store";
import Errors from "@/components/block/Errors";

export default {
  name: 'PostProcessing',
  components: {Errors},
  store: store,
  props: ["selectedRunIds", "plotable"],
  data() {
    return {
      expectationValue: null,
      xAxis: "",
      yAxis: "",
      flags: "",
      selectedParameters: null,

      outputOptions: [
        "Plot Selected",
        "List Inputs",
        "Input Differences"
      ],
      selectedOutputOptions: {
        "Plot Selected": true,
        "List Inputs": true,
        "Input Differences": true
      },
      expecItemList: null, // ['Laser', 'Expec', 'Eig'],
      parameterList: null, //[
      //   {text: "x=0 y=1 -linY", value: {xAxis: 0, yAxis: 1, flags: "-linY"}},
      //   {text: "x=2 y=4 -compare", value: {xAxis: 2, yAxis: 4, flags: "-compare"}},
      //   {text: "x=2 y=3 -linY", value: {xAxis: 2, yAxis: 3, flags: "-linY"}},
      //   {text: "x=5 y=9 -linY", value: {xAxis: 5, yAxis: 9, flags: "-linY"}}
      // ],

      plotImageUrl: null,
      listInputsOutput: null,
      inputDifferenceOutput: null,

      processingPlot: false,
      processingListInputs: false,
      processingInputDifference: false,

      errors: []
    }
  },
  computed: {
    hasPlotSelected() {
      return (this.selectedRunIds && this.selectedRunIds.length === 1) || this.selectedOutputOptions['Plot Selected'];
    },
    hasListInputsSelected() {
      return this.selectedRunIds && this.selectedRunIds.length > 1 && this.selectedOutputOptions['List Inputs'];
    },
    hasInputDifferenceSelected() {
      return this.selectedRunIds && this.selectedRunIds.length > 1 && this.selectedOutputOptions['Input Differences'];
    }
  },
  methods: {
    selectAllOutputOptions() {
      for (let outputOption in this.selectedOutputOptions) {
        this.selectedOutputOptions[outputOption] = true;
      }
    },
    unselectAllOutputOptions() {
      for (let outputOption in this.selectedOutputOptions) {
        this.selectedOutputOptions[outputOption] = false;
      }
    },
    async plotRuns() {
      this.errors = [];
      const promises = [];
      if (this.hasPlotSelected) {
        this.plotImageUrl = null;
        this.processingPlot = true;
        promises.push(PlotService.plotSelectedRuns({
          runIds: this.selectedRunIds,
          expectationValue: this.expectationValue,
          xAxis: this.xAxis,
          yAxis: this.yAxis,
          flags: this.flags
        }).then(({plotImageUrl, output, userGuidance}) => {
          if (plotImageUrl) {
            this.plotImageUrl = plotImageUrl;
          } else {
            const e = new Error();
            e.stack = `[User Guidance]\n${userGuidance}\n[Output]\n${output}`;
            throw e;
          }
          this.processingPlot = false;
        }).catch(e => {
          this.processingPlot = false;
          this.errors.push({
            variant: "danger",
            title: "Error Plotting",
            description: "",
            source: e
          });
        }));
      }

      if (this.hasListInputsSelected) {
        this.listInputsOutput = null;
        this.processingListInputs = true;
        promises.push(PlotService.getRunListInputs({runIds: this.selectedRunIds}).then(({output}) => {
          this.listInputsOutput = output;
          this.processingListInputs = false;
        }).catch(e => {
          this.processingListInputs = false;

          if (e.response && e.response.data) {
            e.stack = `[Output]\n${e.response.data.output}\n[Standard Error]\n${e.response.data.stderr}`;
          }

          this.errors.push({
            variant: "danger",
            title: "Error List Inputs",
            description: "",
            source: e
          });
        }));
      }

      if (this.hasInputDifferenceSelected) {
        this.inputDifferenceOutput = null;
        this.processingInputDifference = true;
        promises.push(PlotService.getRunListInputDifference({runIds: this.selectedRunIds}).then(({output}) => {
          this.inputDifferenceOutput = output;
          this.processingInputDifference = false;
        }).catch(e => {
          this.processingInputDifference = false;

          if (e.response && e.response.data) {
            e.stack = `[Output]\n${e.response.data.output}\n[Standard Error]\n${e.response.data.stderr}`;
          }

          this.errors.push({
            variant: "danger",
            title: "Error Input Differences",
            description: "",
            source: e
          });
        }));
      }

      await Promise.all(promises);
      // this.plotImageUrl = output[0];
      // this.listInputsOutput = output[1];
      // this.inputDifferenceOutput = output[2];

      // So that the expectation value dropdown would have the latest parameters just executed on the top.
      this.expecItemList = await PlotService.getPlotables({runIds: this.selectedRunIds});
    },
    async reset() {
      this.expectationValue = this.plotable;

      this.selectedOutputOptions["Plot Selected"] = true;
      this.selectedOutputOptions["List Inputs"] = true;
      this.selectedOutputOptions["Input Differences"] = true;

      this.plotImageUrl = null;
      this.listInputsOutput = null;
      this.inputDifferenceOutput = null;

      this.parameterList = null;
      this.expecItemList = null;

      this.errors = [];

      await Promise.all([
        PlotService.getPlotables({runIds: this.selectedRunIds}).then(expecItemList => {
          this.expecItemList = expecItemList;
          this.processingExpecItemList = false;

          if (!this.expecItemList || this.expecItemList.length === 0) {
            this.errors.push({
              variant: "warning",
              title: "Invalid Run Selection",
              description: "The selected run(s) do not have any plotables available.",
            });
          }
        }).catch(e => {
          this.processingExpecItemList = false;
          this.errors.push({
            variant: "danger",
            title: "Network Error",
            description: "Unknown error when fetching plotables.",
            source: e
          });
        }),
        PlotService.getPlotParameters().then(parameterList => {
          this.parameterList = parameterList;
          this.processingParameterList = false;
        }).catch(e => {
          this.processingParameterList = false;
          this.errors.push({
            variant: "danger",
            title: "Network Error",
            description: "Unknown error when fetching expectation values.",
            source: e
          });
        })
      ]);
    }
  },
  watch: {
    parameterList() {
      if (this.parameterList && this.parameterList.length > 0) {
        this.selectedParameters = this.parameterList[0].value;
      }
    },
    expecItemList() {
      if (!this.plotable && this.expecItemList && this.expecItemList.length > 0) {
        this.expectationValue = this.expecItemList[0];
      }
    },
    selectedParameters() {
      if (this.selectedParameters) {
        this.xAxis = this.selectedParameters.xAxis;
        this.yAxis = this.selectedParameters.yAxis;
        this.flags = this.selectedParameters.flags;
      }
    },
    selectedRunIds(nw, od) {
      if (JSON.stringify(nw) !== JSON.stringify(od)) {
        this.reset();
      }
    },
    plotable() {
      this.expectationValue = this.plotable;
    },
    expectationValue() {
      this.$emit("plotable-change", this.expectationValue);
    }
  },
  mounted() {
    this.reset();
  }
}
</script>

<style scoped>
.console-output {
  overflow-wrap: anywhere;
  white-space: pre-line;
  background-color: #d1e5cc;
  font-family: "Consolas", monospace;
  font-size: 10px;
}
</style>
