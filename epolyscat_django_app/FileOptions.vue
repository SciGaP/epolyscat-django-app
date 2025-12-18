<template>
    <div class="w-100 d-flex flex-column-reverse justify-content-center" style="padding: 0 20px">
        <!-- TEXT_TABLE This needs to be at the top or else it causes problems and I'm using height instead of v-if due to vue refs (and another weird bug that happens with v-show): 
        https://stackoverflow.com/questions/54355375/vue-js-refs-are-undefined-even-though-this-refs-shows-theyre-there -->
        <LoadingOverlay
            style="width: 100%; position: relative; overflow: hidden"
            :style="{'height': (choice == 'text' && !cantViewContents) ? '400px' : '0'}"
            name="codemirror"
        ><div ref="editor" class="bg-dark" style="height: 400px; width: 100%;"/></LoadingOverlay>
        <!-- NONE -->
        <h2 v-if="choice == 'none'" class="contentBlock" style="opacity: 0.7;">No Preview available</h2>
        <!-- TABLE -->
        <TableView v-if="choice == 'table'" :viewing="runViewType == RUN_VIEW_TYPE.viewing" :selectedTableObject="selectedInfo.tableObject"/>
        <!-- PLOT -->
        <PlotView v-if="choice == 'plot'" :runId="runId"/>
        <!-- UPLOAD -->
        <div v-if="choice == 'upload'" class="w-100 d-flex flex-row">
            <div class="upload-box">
                <b-button variant="primary" v-b-modal:upload-files>Upload from storage</b-button>
                <UserStorage id="upload-files" selectionType="files" :onlyShowEditable="false" @filesSelected="addToInputFilesFromStorage"/>
                <label for="uploadFile" class="btn btn-primary">Upload from computer</label>
                <input style="display: none" type="file" class="btn btn-primary" 
                    id="uploadFile" v-on:change="addToInputFilesFromComputer" multiple>
            </div>
            <div class="upload-files-box">
                <span v-for="filename in uploadedFiles" v-bind:key="filename" > {{ filename }} </span>
            </div>
        </div>
        <!-- Text and Table buttons -->
        <div 
            class="w-100 d-flex flex-row justify-content-between align-items-end" style="margin-bottom: 10px"
        >
            <span>{{ selectedInfo.selected }}</span>
            <div class="d-flex flex-row" v-if="!cantViewContents && choice != 'upload'">
                <b-button 
                    variant="outline-primary" size="sm" :disabled="!isAvailable.text" 
                    @click="choice = 'text'" :pressed="choice == 'text'"
                ><b-icon icon="text-left" /></b-button>
                <b-button 
                    variant="outline-primary" size="sm" :disabled="!isAvailable.table" 
                    @click="choice = 'table'" :pressed="choice == 'table'"
                ><b-icon icon="table" /></b-button>
                <b-button 
                    variant="outline-primary" size="sm" :disabled="!isAvailable.plot"
                    @click="choice = 'plot'" :pressed="choice == 'plot'"
                ><b-icon icon="bar-chart-line" /></b-button>
                <!-- graph-up -->
            </div>
        </div>
        <ReplaceFileModal ref="replaceFileModal"/>
    </div>
</template>

<script>
import UserStorage from '../overlay/UserStorage.vue';
import TableView from './TableView.vue';
import LoadingOverlay from '../overlay/LoadingOverlay.vue';
import { eventBus } from '@/event-bus';
import store from '@/store';
import CodeMirror from "codemirror";
import "codemirror/lib/codemirror.css";
import "codemirror/theme/material.css";
import PlotView from './PlotView.vue';
import ReplaceFileModal from '../blocks/ReplaceFileModal.vue';

export default {
    props: ["runViewType", "runId"],
    store,
    data() {
        return {
            editor: null,
            contents: null,
            fileIsLoaded: true,
            choice: "text",
        };
    },
    components: { UserStorage, TableView, LoadingOverlay, PlotView, ReplaceFileModal },
    created() {
        // TODO: scrap this and do something better, perhaps move to seperate file?
        this.VIEW_TYPE = { 
            none: "none", 
            upload: "upload", 
            text: "text", 
            table: "table",
            plot: "plot"
        };

        this.RUN_VIEW_TYPE = {
            creating: "creating",
            editing: "editing",
            viewing: "viewing"
        }
    },
    computed: {
        selectedInfo() {
            return this.$store.getters["input/getSelectedInfo"]();
        },
        properFilename() {
            return this.$store.getters["input/getProperFilename"];
        },
        uploadedFiles() {
            const { inputFiles } = this.$store.getters["input/getInputs"]();
            const path = this.$store.getters["input/getPath"];

            if (path.length == 4)
                return inputFiles[path[2]].files.map(file => file.name) 
            
            return Object
                .keys(inputFiles)
                .filter(key => inputFiles[key].files.length > 0)
        },
        path() {
            return this.$store.getters["input/getPath"];
        },
        cantViewContents() {
            return this.contents == null;
        },
        isAvailable() {
            return {
                'plot': this.selectedInfo.plotObject != null,
                'table': this.selectedInfo.tableObject != null,
                'text': !this.cantViewContents,
                'upload': this.selectedInfo.selected == 'Upload',
                'none': true,
            };
        },
    },
    methods: {
        addToInputFilesFromComputer(event) {
            let files = Array.from(event.target.files);
            files.forEach(file => file.isFromComputer = true);

            this.addToInputFiles(files);
        },
        addToInputFilesFromStorage(selected) {
            selected.forEach(file => file.isFromComputer = false);

            this.addToInputFiles(selected);
        },
        async addToInputFiles(files) {
            // let preExistingFiles = this.$store.getters["input/getFilesInCategory"](this.path[this.path.length-2])
            //     .filter(file => !file.deleted);
            let prevAnswer = null;

            // if (!this.$store.getters["input/multipleFilesCanBeUploaded"]() && (preExistingFiles.length > 0 || files.length > 1)) {
            //     let _ = await this.$bvModal.msgBoxOk(
            //         `"${this.path[this.path.length-2]}" can only have one file uploaded to it`, {
            //             centered: true
            //         }
            //     );

            //     if (preExistingFiles.length > 0)
            //         return;
            //     else
            //         files = files.slice(0, 1);
            // }

            for (let file of files) {
                if (this.$store.getters["input/checkIfDuplicate"](file.name)) {
                    const answer = 
                        (prevAnswer == "none" || prevAnswer == "all") ? prevAnswer 
                        : await this.$refs.replaceFileModal.ask(file.name);

                    prevAnswer = answer;

                    file.replaceCurrent = (answer == "all" || answer == "yes");
                }

                this.$store.dispatch("input/addFiles", { files: [file] });
            }
        },
        updateParameters(parameters) {
            this.$store.commit("input/SET_INPUTS", { parameters });
        },
        loadFileContents(filename) {
            this.$store.commit("loading/START", { key: "codemirror", message: "Loading File" });
            // Prevents selecting files in rapid succession causing 
            // the editor to switch between contents before settling 
            // on the correct contents

            this.fileIsLoaded = false;
            this.contents = "";

            try {
                this.$store.getters["input/getContentsOfFile"]()
                    .then(contents => {
                        if (filename == this.selectedInfo.selected && !this.fileIsLoaded) {
                            this.contents = contents
                            this.fileIsLoaded = true;

                            if (contents != null) {
                                this.editor.setValue(contents);
                            }

                            this.$store.commit("loading/STOP", { key: "codemirror", message: "Loading File" });
                        }
                    });
            } catch (error) { 
                eventBus.$emit("error", { name: `Error while trying to get contents of "${this.selected}"`, error }); 
            }
        },
        updateChoice() {
            const { selected } = this.selectedInfo;

            if (this.isAvailable[this.choice] === false || this.choice == 'none') 
                this.choice = Object.entries(this.isAvailable).find(([_, isAvailable]) => isAvailable)[0];
        }
    },
    watch: {
        async choice() {
            const { selected, isFile } = this.selectedInfo;

            if (this.choice == "text" && isFile) {
                this.loadFileContents(selected);
            }
        },
        path() {
            const { selected, isFile } = this.selectedInfo;

            if (isFile) 
                this.loadFileContents(selected);
            else
                this.contents = null;

            this.updateChoice();
        },
        isAvailable() {
            this.updateChoice();
        },
        runViewType(runViewType) {
            this.editor.options.readOnly = (runViewType == this.RUN_VIEW_TYPE.viewing);
        },
        contents() {
            this.updateChoice();
        }
    },
    mounted() {
        this.editor = CodeMirror(this.$refs.editor, {
            theme: "material",
            mode: "text/plain",
            lineNumbers: true,
            lineWrapping: false,
            scrollbarStyle: "native",
            extraKeys: {"Ctrl-Space": "autocomplete"},
            value: ""
        });

        this.editor.setSize("100%", "100%");
        this.editor.on("change", () => {
            window.setTimeout(() => this.editor.refresh(), 1);

            if (this.contents != null) {
                this.$store.dispatch("input/setContents", { 
                    contents: this.editor.getValue() 
                });
            }
        });
    }
}
</script>

<style scoped>

    .contentBlock {
        display: flex;
        justify-content: center;
        align-items: center;
        border: 1px solid #ddd; 
        border-radius: 8px;
        width: 100%; 
        height: 400px;
    }

    /* Upload box */

    .upload-box {
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        margin: 0;
        border: 4px dashed #ddd;
        border-radius: 8px;
        height: 400px;
        flex-grow: 1;
    }
    
    .upload-files-box {
        display: flex;
        flex-direction: row;
        flex-wrap: wrap;
        align-content: baseline;
        margin: 0 10px;
        border-radius: 8px;
        padding: 5px 10px;
        height: 400px;
        max-width: 40%;
        background: #eee;
        flex-grow: 1;
    }

    .upload-box > * {
        margin: 10px 20px;
    }

    .upload-files-box > * {
        margin: 10px;
        height: fit-content;
    }
</style>