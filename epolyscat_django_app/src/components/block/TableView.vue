<template>
    <LoadingOverlay :name="loaderName">
        <b-card no-body>
            <b-tabs pills card content-class="w-100 overflow-auto" nav-class="pages-nav-pills">
                <b-tab :title="page.name" v-for="page in tableObject.pages" :key="page" no-body>
                    <div class="tableView">
                        <template v-if="page.type == 'table'">
                            <DataInput 
                                v-for="item in filterItems(page.data)" :key="{dataKey, item}" 
                                :item="item" style="min-width: 120px; margin: 20px; max-width: 400px;"
                                :includeName="true" :viewing="viewing"
                            />
                        </template>
                        <template v-else-if="page.type == 'list'">
                            <b-table-simple sticky-header="400px" responsive style="height: 390px; margin-bottom: 0;" ref="simpleTable">
                                <b-thead>
                                    <b-tr>
                                        <b-th 
                                            v-for="column in page.columns" 
                                            :key="column" class="text-nowrap pl-3"
                                            style="background: #d4e0ea; color: #226597;"
                                        >
                                            {{ column.name }}
                                        </b-th>
                                        <b-th v-if="!viewing" style="background: #d4e0ea;"/>
                                    </b-tr>
                                </b-thead>
                                <b-tbody>
                                    <b-tr v-for="(row, i) in page.data" :key="row" style="border-top: 1px solid gray">
                                        <b-td
                                            v-for="(column, j) in page.columns" 
                                            :key="column" class="text-nowrap"
                                        >
                                            <div class="w-100 h-100 d-flex">
                                                <DataInput 
                                                    :item="{...column, ...row.cells[j]}" :key="dataKey"
                                                    :includeName="false" :viewing="viewing" v-if="j in row.cells"
                                                />
                                                <div v-else-if="!viewing && (j == 0 || (j - 1) in row.cells)">
                                                    <b-button v-if="'removeCell' in row && !(i == 0 && j == 0) " variant="link" size="sm" @click="removeCell(row)">
                                                        Delete cell
                                                    </b-button>
                                                    <b-button v-if="'addCell' in row" variant="link" size="sm" @click="addCell(row)">
                                                        Add cell
                                                    </b-button>
                                                </div>
                                            </div>
                                        </b-td>
                                        <b-td v-if="!viewing" style="vertical-align: middle">
                                            <b-button variant="link" size="sm" @click="deleteRow(row)">
                                                <b-icon icon="trash"/>
                                            </b-button>
                                        </b-td>
                                    </b-tr>
                                </b-tbody>
                                <b-tfoot>
                                    <b-tr v-if="!viewing && 'addRow' in page">
                                        <b-td style="vertical-align: middle">
                                            <b-button variant="link" size="sm" @click="addRow(page)">
                                                Add a new row
                                            </b-button>
                                        </b-td>
                                        <b-td :colspan="page.columns.length"/>
                                    </b-tr>
                                </b-tfoot>
                            </b-table-simple>
                        </template>
                    </div>
                </b-tab>
            </b-tabs>
        </b-card>
    </LoadingOverlay>
</template>

<script>
import DataInput from './DataInput.vue';
import store from '@/store';
import LoadingOverlay from '../overlay/LoadingOverlay.vue';
import { eventBus } from '@/event-bus';

export default {
    components: { DataInput, LoadingOverlay },
    store,
    props: {
        viewing: Boolean,
        showParameters: Boolean,
        selectedTableObject: Object,
        key: {
            type: Number,
            required: false
        }
    },
    data() {
        return {
            tableObject: { pages: [] },
            dataKey: 0
        };
    },
    computed: {
        loaderName() {
            return `table-${this._uid}`;
        }
    },
    methods: {
        async loadTableObject() {
            // console.log("selected table obj", this.selectedTableObject);
            this.$store.commit("loading/START", { key: this.loaderName, message: "Loading Table" });
            
            const scrollTop = ("simpleTable" in this.$refs) ? this.$refs.simpleTable[0].$el.scrollTop : null;

            try {
                if (this.selectedTableObject == null)
                    throw Error(`Could not get the table object for "${selected}"`);

                this.tableObject = await this.getTableObjectData(this.selectedTableObject);

                this.dataKey += 1;

                if (scrollTop)
                    this.$nextTick(function () { this.$refs.simpleTable[0].$el.scrollTo(0, scrollTop); });

                this.$store.commit("loading/STOP", { key: this.loaderName, message: "Loading Table" });
            } catch (error) {
                eventBus.$emit("error", { name: "Error while trying to load table", error });
            }
        },
        async getTableObjectData(tableObject) {
            let newTableObject = (Array.isArray(tableObject)) ? [] : {};

            if ("get" in tableObject) {
                const data = await tableObject.get();

                if (typeof data == "object" && data != null)
                    newTableObject.data = await this.getTableObjectData(data);
                else
                    newTableObject.data = data;
            }

            for (const key in tableObject) {
                if (typeof tableObject[key] == "object" && tableObject[key] != null)
                    newTableObject[key] = await this.getTableObjectData(tableObject[key]);
                else
                    newTableObject[key] = tableObject[key];
            }

            return newTableObject;
        },
        filterItems(items) {
            return items.filter(item => 'data' in item);
        },
        async addRow(page) {
            await page.addRow(); 
            await this.loadTableObject();
        },
        async removeCell(row) {
            await row.removeCell();
            await this.loadTableObject();
        },
        async deleteRow(row) {
            await row.delete();
            await this.loadTableObject();
        },
        async addCell(row) {
            await row.addCell();
            await this.loadTableObject();
        },
    },
    watch: {
        selectedTableObject() {
            this.loadTableObject();
        }
    },
    mounted() {
        // console.log("MOUNTED AND GOING", this);
        this.loadTableObject();
    }
}
</script>

<style lang="scss">
@import "../../styles";

.tableView {
    display: flex;
    flex-direction: row;
    flex-wrap: wrap;
    align-content: flex-start;
    overflow: scroll;
    min-width: 100%; 
    max-height: 450px;
    height: fit-content;
}
    
//     ul.pages-nav-pills li.nav-item a.nav-link,
// ul.sections-nav-pills li.nav-item a.nav-link {
//     border-radius: 20px;
//     padding: 5px 30px;
//     margin-right: 10px;
//     margin-bottom: 10px;
//     background-color: $light;
// }

// ul.pages-nav-pills li.nav-item a.nav-link.active {
//     background-color: $primary;
//     color: $light;
// }

// ul.pages-nav-pills li.nav-item a.nav-link {
//     background-color: $light;
//     color: $primary;
// }

// ul.sections-nav-pills li.nav-item a.nav-link.active {
//     background-color: $info;
//     color: $light;
// }

// ul.sections-nav-pills li.nav-item a.nav-link {
//     background-color: $light;
//     color: $info;
// }
</style>