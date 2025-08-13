<template>
  <div class="w-100 p-3">
    <Errors :errors="errors"/>
    <div class="w-100 pb-2 d-flex flex-row">
      <b-form-select v-if="filenames && filenames.length > 0" v-model="filename" :options="filenames"
                     style="max-width: 100px;" size="sm"/>
      <b-button v-if="allowRefresh" size="sm" variant="link" title="Refresh viewable content" v-b-tooltip.hover
                v-on:click="refreshData" :disabled="processing">
        <template v-if="processing">
          <b-icon icon="arrow-clockwise" animation="spin"></b-icon>
          Refreshing ...
        </template>
        <template v-else>
          <b-icon icon="arrow-clockwise"></b-icon>
          Refresh
        </template>
      </b-button>
      <div class="text-right flex-fill pr-3">Switch input view:</div>
      <b-button size="sm" :variant="inputType === 'text' ? 'primary': 'outline-primary'" :disabled="!run"
                v-on:click="inputType='text'" v-b-tooltip.hover title="Text input">
        <b-icon icon="file-text"/>
      </b-button>
      <b-button size="sm" :variant="inputType === 'table' ? 'primary': 'outline-primary'"
                :disabled="!run || !run.inputTable"
                v-on:click="inputType='table'" v-b-tooltip.hover title="Table input">
        <b-icon icon="table"/>
      </b-button>
    </div>
    <b-overlay :show="inputType === 'table' && runId && (!run || !run.inputTable)" rounded="sm">
      <RunInputTable v-if="run" :class="{visible: inputType === 'table', invisible: inputType !== 'table'}"
                     class="w-100" :inputTable="run.inputTable" v-on:change="onInputTableValueChange"
                     :read-only="readOnly"/>
    </b-overlay>
    <b-overlay :show="inputType === 'text' && runId && (!run || !editorValue)" rounded="sm">
      <div :class="{visible: inputType === 'text', invisible: inputType !== 'text'}" class="w-100 bg-dark" ref="editor"
           style="min-height: 300px;">
      </div>
    </b-overlay>
  </div>
</template>

<script>
import CodeMirror from "codemirror";
import "codemirror/lib/codemirror.css";
import "codemirror/theme/abcdef.css";
import store from "@/store";
import RunInputTable from "@/components/block/run-input/RunInputTable";
import Errors from "@/components/block/Errors";

export default {
  name: 'RunViewableEditor',
  components: {Errors, RunInputTable},
  store: store,
  props: {
    "runId": {
      default: null
    },
    "run": {
      default: null
    },
    "filename": {
      default: "inpc" // Default input file.
    },
    "filenames": {
      default: null
    },
    "readOnly": {
      default: true
    },
    "allowRefresh": {
      default: false
    },
  },
  data() {
    return {
      inputType: 'text',
      editor: null,
      value: null,
      selectedFilename: null,

      processing: false,
      errors: []
    };
  },
  computed: {
    editorValue() {
      let _editorValue;
      if (this.run) {
        _editorValue = this.$store.getters["run/getViewableContent"]({
          runId: this.run.runId,
          filename: this.filename,
          inpcDownloadUrl: this.run.inpcDownloadUrl
        });
      }

      if (_editorValue) {
        return _editorValue;
      } else {
        return "";
      }
    }
  },
  methods: {
    async refreshData() {
      if (this.run) {
        this.errors = [];
        this.processing = true;

        this.editor.setValue(this.editorValue);

        try {
          await this.$store.dispatch("run/fetchViewableContent", {
            runId: this.run.runId,
            filename: this.filename,
            inpcDownloadUrl: this.run.inpcDownloadUrl
          });
        } catch (e) {
          this.errors.push({
            variant: "danger",
            title: "Error List Inputs",
            description: "Unknown error when creating inputs list.",
            source: e
          });
        }

        this.processing = false;
      }
    },
    onInputTableValueChange(value) {
      this.value = value;
    }
  },
  watch: {
    editorValue() {
      if (this.editorValue) {
        this.editor.setValue(this.editorValue);
      }
    },
    run() {
      this.refreshData();
    },
    filename() {
      this.selectedFilename = this.filename;
      this.refreshData();
    },
    value() {
      this.$emit("change", this.inputType, this.value);
    },
    inputType() {
      if (this.inputType === "text") {
        this.value = this.editor.getValue();
      } else if (this.inputType === "table" && this.run) {
        this.value = this.run.inputTable;
      }
    }
  },
  async mounted() {
    this.selectedFilename = this.filename;

    this.editor = CodeMirror(this.$refs.editor, {
      theme: "abcdef",
      mode: "text/plain",
      readOnly: this.readOnly,
      lineNumbers: true,
      lineWrapping: false,
      scrollbarStyle: "native",
      extraKeys: {"Ctrl-Space": "autocomplete"},
      value: "",
    });
    this.editor.setSize("100%", "100%");

    this.editor.on("change", () => this.value = this.editor.getValue());

    this.refreshData();
  }
};
</script>

<style scoped>
.visible {
  visibility: unset;
}

.invisible {
  visibility: hidden;
  position: fixed;
  top: -10000px;
}
</style>
