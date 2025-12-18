<template>
    <div class="h-100 d-flex flex-column" style="width: 275px">
        <router-link :to="`/`" v-slot="{ href, navigate, isExactActive }" class="link">
            <b-link variant="link" class="d-flex flex-row text-center" :class="{active: isExactActive}" :href="href" @click="navigate">
                <b-icon font-scale="1.2" icon="house-door"/>
                <span>Home</span>
            </b-link>
        </router-link>
        <div class="d-flex flex-row">
            <router-link :to="`/runs`" v-slot="{ href, navigate, isExactActive }" class="link">
                <b-link variant="link" class="d-flex flex-row text-center" :class="{active: isExactActive}" :href="href" @click="navigate">
                    <b-icon font-scale="1.2" icon="collection"/>
                    <span>Runs</span>   
                </b-link>
            </router-link>
            <b-button 
                variant="link" size="sm" aria-controls="moreInfo" class="ml-2"
                v-b-toggle.show-runs
            >
                <b-icon :icon="(listRuns) ? 'chevron-down' : 'chevron-right'"></b-icon>
            </b-button>
        </div>
        <b-collapse id="show-runs" v-model="listRuns" accordion="show-more" role="tabpanel">
            <div style="border-left: 1px solid black; margin-left: 40px;">
                <LoadingOverlay name="runs" class="d-flex flex-column">
                    <router-link style="margin-bottom: 0;" v-for="run in runs" v-bind:key="run" :to="`/runs/${run.id}`" v-slot="{ href, navigate, isExactActive }" class="link">
                        <b-link variant="link" class="d-flex flex-row text-center" :class="{active: isExactActive}" :href="href" @click="navigate">
                            <span class="ellipses">{{ run.name }}</span>
                        </b-link>
                    </router-link>
                    <div v-if="runs.length == 0" style="color: gray; margin: 20px">No Runs</div>
                    <router-link v-else :to="`/runs`" v-slot="{ href, navigate, isExactActive }" class="link">
                        <b-link variant="link" class="d-flex flex-row text-center" :class="{active: isExactActive}" :href="href" @click="navigate">
                            <span style="margin: 0 5px">All Runs</span>
                        </b-link>
                    </router-link>
                </LoadingOverlay>
            </div>
        </b-collapse>
        <div class="d-flex flex-row">
            <router-link :to="`/views`" v-slot="{ href, navigate, isExactActive }" class="link">
                <b-link variant="link" class="d-flex flex-row text-center" :class="{active: isExactActive}" :href="href" @click="navigate">
                    <b-icon font-scale="1.2" icon="eye"/>
                    <span>Views</span>
                </b-link>
            </router-link>
            <b-button 
                variant="link" size="sm" aria-controls="moreInfo" class="ml-2"
                v-b-toggle.show-views
            >
                <b-icon :icon="(listViews) ? 'chevron-down' : 'chevron-right'"></b-icon>
            </b-button>
        </div>
        <b-collapse id="show-views" v-model="listViews" accordion="show-more" role="tabpanel">
            <div style="border-left: 1px solid black; margin-left: 40px">
                <LoadingOverlay name="views" class="d-flex flex-column">
                    <router-link style="margin-bottom: 0;" v-for="view in views" v-bind:key="view" :to="`/views/${view.id}`" v-slot="{ href, navigate, isExactActive }" class="link">
                        <b-link variant="link" class="d-flex flex-row text-center" :class="{active: isExactActive}" :href="href" @click="navigate">
                            <span class="ellipses">{{ view.name }}</span>
                        </b-link>
                    </router-link>
                    <div v-if="views.length == 0" style="color: gray; margin: 20px">No Views</div>
                    <router-link v-else :to="`/views`" v-slot="{ href, navigate, isExactActive }" class="link">
                        <b-link variant="link" class="d-flex flex-row text-center" :class="{active: isExactActive}" :href="href" @click="navigate">
                            <span style="margin: 0 5px">All Views</span>
                        </b-link>
                    </router-link>
                </LoadingOverlay>
            </div>
        </b-collapse>
        <router-link :to="`/tutorials`" v-slot="{ href, navigate, isExactActive }" class="link">
            <b-link variant="link" class="d-flex flex-row text-center" :class="{active: isExactActive}" :href="href" @click="navigate">
                <b-icon font-scale="1.2" icon="collection"/>
                <span>Tutorials</span>   
            </b-link>
        </router-link>
    </div>
</template>

<script>
import store from '@/store';
import LoadingOverlay from '../overlay/LoadingOverlay.vue';

export default {
    store,
    components: { LoadingOverlay },
    data() {
        return {
            listRuns: true,
            listViews: false
        }
    },
    computed: {
        views() {
            let views = this.$store.getters["view/getViews"](); 

            views.sort((view1, view2) => (new Date(view2.updated)).getTime() - (new Date(view1.updated)).getTime())
            
            return views.slice(0, 4);
        },
        runs() {
            let runs = this.$store.getters["run/getRuns"](); 

            runs.sort((run1, run2) => (new Date(run2.updated)).getTime() - (new Date(run1.updated)).getTime())
            
            return runs.slice(0, 4);
        },
    },
    methods: {
    },
    async mounted() {
        // this.$store.commit("loading/START", { key: "runs", message: "Fetching Runs" });
        // this.$store.commit("loading/START", { key: "views", message: "Fetching Views" });

        // try {
        //     await this.$store.dispatch("run/fetchRuns", {});

        //     this.$store.commit("loading/STOP", { key: "runs", message: "Fetching Runs" });
        // } catch (error) { 
        //     eventBus.$emit("error", { name: `Error while trying to fetch the runs`, error }); 
        // }

        // try {
        //     await this.$store.dispatch("run/fetchRuns", {});

        //     this.$store.commit("loading/STOP", { key: "runs", message: "Fetching Views" });
        // } catch (error) { 
        //     eventBus.$emit("error", { name: `Error while trying to fetch the views`, error }); 
        // }
    },
}
</script>

<style scoped>
    .link > * {
        margin: 0 5px;
    }

    .link {
        margin: 20px;
        line-height: normal;
    }

    .ellipses {
        text-overflow: ellipsis;
        white-space: nowrap;
        overflow: hidden;
    }
</style>
