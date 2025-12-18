<template>
    <div class="d-flex flex-row justify-content-center flex-wrap">
        <ButtonsGroup 
            v-if="!viewing" label="Run Type/App Module" @setPath="setPathAt"
            :buttonsList="getButtonsList().slice(0, 2)" :offset="0" class="mr-5 mb-2"
        />
        <ButtonsGroup 
            v-if="!showOutputs" label="Input Files" @setPath="setPathAt" 
            :buttonsList="getButtonsList().slice(2)" :offset="2" class="mb-2"
        />
        <ButtonsGroup 
            v-if="showOutputs" label="Files" @setPath="setPathAt"
            :buttonsList="getButtonsList().slice(2)" :offset="2" class="mb-2"
        />
    </div>
</template>

<script>
import store from "@/store";
import ButtonsGroup from "./ButtonsGroup.vue";

export default {
    store,
    props: ["viewing"],
    components: { ButtonsGroup },
    computed: {
        displayedButtonsList() {
            return this.getButtonsList().slice(this.offset);
        },
        offset() {
            return (this.viewing) ? 2 : 0
        },
        path() {
            return this.$store.getters["input/getPath"];
        },
        showOutputs() {
            return Object.keys(this.$store.getters["input/getInputs"]().outputFiles).length > 0
        }
    },
    methods: {
        async loadPathLabels() {
            await this.$store.dispatch("input/fetchPathLabels", {});

            this.selectDefaultPathFrom(0);
        },
        setPathAt({index, key}) {
            let path = this.path.slice();
            path[index] = key;

            if (this.viewing || index >= 2) this.$emit("clickedPath");

            this.$store.commit("input/SET_PATH", { path });
            this.selectDefaultPathFrom(index + 1) 
        },
        selectDefaultPathFrom(index) {
            let path = this.path.slice(0, index);
            let buttonsList = this.getButtonsList({ customPath: path });

            while (buttonsList.length > path.length) {
                if (buttonsList.length == 4 && buttonsList[3].length > 1 && !this.viewing)
                    path.push(buttonsList[3][1].key);
                else if (buttonsList[buttonsList.length - 1].length > 0)
                    path.push(buttonsList[buttonsList.length - 1][0].key);
                else 
                    break;

                buttonsList = this.getButtonsList({ customPath: path });
            }

            this.$store.commit("input/SET_PATH", { path });
        },
        getButtonsList({ customPath = null } = {}) {
            const {
                pathLabels,
                inputFiles, 
                parameters,
                outputFiles
            } = this.$store.getters["input/getInputs"]({ customPath });

            const files = (Object.keys(outputFiles).length > 0) ? outputFiles : inputFiles;

            const path = (customPath != null) ? customPath : this.path;

            let buttonsList = [[
                { key: "MODULE" },
                { key: "UTILITY" },
                { key: "WORKFLOW" },
            ]];

            if (path.length >= 1)
                buttonsList[1] = pathLabels[1].map(key => ({ key }));

            if (path.length >= 2) {
                let inputsButtons = [];
                
                if (!this.viewing)
                    inputsButtons.push({ key: "Upload" });

                inputsButtons = inputsButtons.concat(Object.entries(files).map(([
                    fileGroupName,
                    file
                ]) => ({
                    key: fileGroupName,
                    checked: file.isValid(inputFiles)
                })));

                console.log("Parameters: ", parameters);

                // if (parameters.length > 0)
                //     inputsButtons.push({ 
                //         key: "Parameters", 
                //         checked: parameters.every(parameter => parameter.issues(parameter.value).length == 0)
                //     });

                if (this.viewing)
                    inputsButtons = inputsButtons.map(({ key }) => ({ key }));
                
                buttonsList.push(inputsButtons);
            }

            // console.log(path.length >= 3, path[2] in files);

            if (path.length >= 3 && path[2] in files) {
                let fileButtons = [];

                if (!this.viewing)
                    fileButtons.push({ key: "Upload" });

                if (Array.isArray(files[path[2]].files))
                    fileButtons = fileButtons.concat(files[path[2]].files.map(file => ({
                        key: file.name,
                        delete: () => {
                            this.$store.commit("input/REMOVE_FILE", { 
                                inputFileName: path[2],
                                filename: file.name
                            });
            
                            this.setPathAt({index: 3, key: "Upload"});
                        }
                    })));
                else 
                    fileButtons = fileButtons.concat(Object
                        .keys(files[path[2]].files)
                        .map(categoryName => ({
                            key: categoryName
                        }))
                    );

                if (this.viewing)
                    fileButtons = fileButtons.map(({ key }) => ({ key }));

                buttonsList.push(fileButtons)
            }

            if (path.length >= 4 && path[2] in files && path[3] in files[path[2]].files) {
                let fileButtons = [];
                
                fileButtons = fileButtons.concat(files[path[2]].files[path[3]].files.map(file => ({
                    key: file.name,
                    delete: () => {
                        this.$store.commit("input/REMOVE_FILE", { 
                            inputFileName: path[2],
                            filename: pathOption
                        });
        
                        this.setPathAt({index: 3, key: "Upload"});
                    }
                })));

                if (this.viewing)
                    fileButtons = fileButtons.map(({ key }) => ({ key }));

                buttonsList.push(fileButtons)
            }

            for (const buttons of buttonsList) {
                buttons.perColumn = Math.max(Math.ceil(buttons.length / 2), 5)
            }

            // console.log("buttonsList", buttonsList);

            return buttonsList;
        },
    },
    watch: {
        viewing(viewing) {
            if (viewing) {
                this.selectDefaultPathFrom(2);
            }
        }
    },
    mounted() {
        this.loadPathLabels();
    },
}
</script>
