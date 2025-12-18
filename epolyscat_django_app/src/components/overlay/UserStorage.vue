<template>
    <b-modal size="xl" :id="id" scrollable :title="title" @show="fetchData" body-class="p-0">
        <b-breadcrumb class="my-3 ml-4">
            <b-breadcrumb-item 
                v-for="(directory, i) in path" v-bind:key="directory" 
                @click="breadcrumbClick(i)" 
                :active="i == path.length - 1"
            >{{ directory }}</b-breadcrumb-item>
        </b-breadcrumb>
        <LoadingOverlay name="directory" #default="{  }" style="min-height: 100px; height: calc(100% - 50px)">
            <ListView 
                :items="items" :columns="[['name', 'Name'], ['modifiedTime', 'Last Modified'], ['size', 'Size']]" 
                :canSelectMultiple="selectionType != 'folder'" identifier="name"
                style="height: 100%" @updateSelected="updateSelected" 
                :sorters="[
                    (item1, item2) => item1.name.localeCompare(item2.name), 
                    (item1, item2) => (new Date(item2.createdTime)).getTime() - (new Date(item1.createdTime)).getTime(),
                    (item1, item2) => item1.actualSize - item2.actualSize
                ]"
            >
                <template v-slot:name="{ item }" >
                    <b-link 
                        v-if="item.isDirectory" @click="path.push(item.name)" 
                        style="line-height: 20px; padding-right: 20px; word-break: break-all;"
                    >
                        <b v-if="path.length == 1 && item.name == 'shared'">shared </b>
                        <template v-else>{{ item.name }} </template>
                        <b-icon icon="folder" />
                    </b-link> 
                    <span v-else>{{ item.name }}</span>
                </template>
            </ListView>
        </LoadingOverlay>
        <template #modal-footer="{ ok, cancel }">
            <span class="cutoffText" style="width: 300px">{{ selectedPreview }}</span>
            <div>
                <b-button 
                    variant="outline-primary" :disabled="openDisabled" 
                    @click="disableOpen = true; ok(); emitSelectedFiles(); disableOpen = false;"
                >Open</b-button>
                <b-button variant="outline-dark" @click="cancel()">Cancel</b-button>
            </div>
        </template>
    </b-modal>
</template>

<script>
    import { DirectoryService } from "@/service/epolyscat-service";
    import ListView from "./ListView.vue";
    import LoadingOverlay from "./LoadingOverlay.vue";
    import { eventBus } from "@/event-bus";

    export default {
        props: {
            selectionType: String,
            id: String,
            onlyShowEditable: {
                type: Boolean,
                default: false
            }
        },
        components: { ListView, LoadingOverlay },
        data() {
            return {
                path: ["~"],
                selected: [],
                currentDirectory: {
                    path: [],
                    items: []
                },
                openDisabled: false
            }
        },
        computed: {
            items() {
                let { path, items } = this.currentDirectory;

                // deep copy items
                items = JSON.parse(JSON.stringify(items));

                if (path.join("/") != this.path.join("/"))
                    return [];

                if (this.onlyShowEditable)
                    items = items.filter(item => item.userHasWriteAccess);

                for (let item of items) {
                    item.createdTime = (new Date(item.createdTime)).toLocaleDateString();
                    item.modifiedTime = (new Date(item.createdTime)).toLocaleDateString();
                    
                    item.actualSize = item.size;
                    const orderOfMagnitude = (item.size > 0) ? Math.floor(Math.log(item.size) / Math.log(1000)) : 0;
                    const unit = (orderOfMagnitude == 0) ? "B" : (orderOfMagnitude == 1) ? "KB" : (orderOfMagnitude == 2) ? "MB" : "GB";
                    const size = item.size / (1000)**orderOfMagnitude;
                    item.size = `${Math.floor(100 * size) / 100} ${unit}`;
                    
                    item.isSelectable = (this.selectionType == "folder") ? item.isDirectory : !item.isDirectory;
                }

                if (path.length == 1)
                    items.sort((item1, item2) => (item1.name == "shared") ? -1 : (item2.name == "shared") ? 1 : 0);

                return items;
            },
            isLoaded() {
                return this.currentDirectory.path.join("/") == this.path.join("/");
            },
            selectedPreview() {
                let selected = this.selected;
                
                if (this.selectionType == "folder" && selected.length == 0) {
                    selected = [{
                       name: this.path[this.path.length - 1]
                    }];
                }
                
                return selected.map(file => file.name).join(", ");
            },
            title() {
                return ((this.selectionType == "folder") ? "Select a folder" : "Select all files");
            }
        }, 
        methods: {
            async fetchDirectory(path) {
                this.$store.commit("loading/START", { key: "directory", message: "Loading Directory" });

                try {
                    let {
                        directories, 
                        files
                    } = await DirectoryService.fetchDirectories(path.join("/"));

                    directories = directories.map(directory => ({ isDirectory: true, ...directory }));
                    files = files.map(file => ({ isDirectory: false, ...file }));
                    if (path.join("/") == this.path.join("/")) {
                        this.currentDirectory = {
                            path: path.join("/").split("/"), // copies the path so there is no reference
                            items: [...directories, ...files]
                        }
                    }
                } catch (error) {
                    eventBus.$emit("error", { 
                        name: `Error while trying to load directory: "${path.join("/")}"`, 
                        error 
                    })

                    if (path.join("/") == this.path.join("/")) {
                        this.currentDirectory = {
                            path: path.join("/").split("/"),
                            items: []
                        }
                    }
                }

                this.$store.commit("loading/STOP", { key: "directory", message: "Loading Directory" });
            },
            updateSelected(selected) {
                this.selected = selected;
            },
            emitSelectedFiles() {
                if (this.selectionType == "folder" && this.selected.length == 0) {
                    this.selected = [{
                        path: this.path.join("/")
                    }];
                }

                this.$emit("filesSelected", this.selected);
            },
            breadcrumbClick(i) {
                if (i < this.path.length - 1) {
                    this.path = this.path.slice(0, i + 1);
                    this.selected = [];
                }
            },
            disableOpen() {
                this.openDisabled = true;
            }
        },
        watch: {
            path(path) {
                this.fetchDirectory(path);
                this.selected = [];
            }
        },
        mounted() {
            // console.log("UserStorage mounted");
            this.fetchDirectory(this.path);
        }
    }
</script>
