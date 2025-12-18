<template>
  <div id="app" class="flex-fill h-100 overflow-auto">
    <notifications-display/>
    <div class=" w-100 h-100 d-flex">
      <div class="overflow-auto app_left_nav" style="max-width: 250px;min-width: 250px;box-shadow: 1px 1px 3px 1px #d8d8d8;">
        <AppLeftNav/>
      </div>
      <div class="flex-fill overflow-auto" >
        <router-view class="w-100" :key="$route.fullPath"/>
      </div>
      <!--      <div class="overflow-auto">-->
      <!--        <AppRightNav/>-->
      <!--      </div>-->
    </div>
  </div>
</template>

<script>
import "./styles.scss";
import AppLeftNav from "@/components/blocks/AppLeftNav";
import { eventBus } from "./event-bus";
import ErrorToast from './components/overlay/ErrorToast.vue';
import store from "./store";

const {NotificationsDisplay} = window.CommonUI || {};

export default {
  name: 'App',
  components: {AppLeftNav, NotificationsDisplay},
  store,
  methods: {
        async refreshRuns() {
            this.$store.commit("loading/START", { key: "runs", message: "Fetching Runs" });

            try {
                await this.$store.dispatch("run/fetchRuns", {});

                this.$store.commit("loading/STOP", { key: "runs", message: "Fetching Runs" });
            } catch (error) {
                eventBus.$emit("error", { name: `Error while trying to fetch the runs`, error });
            }
        },
        async refreshViews() {
            this.$store.commit("loading/START", { key: "views", message: "Fetching Views" });

            try {
                await this.$store.dispatch("view/fetchViews", {});

                this.$store.commit("loading/STOP", { key: "views", message: "Fetching Views" });
            } catch (error) {
                eventBus.$emit("error", { name: `Error while trying to fetch the views`, error });
            }
        }
  },

  mounted() {
        eventBus.$on("error", ({ name, error }) => {
            this.$toast.error({
                component: ErrorToast,
                props: { name, error }
            }, {
                position: "top-left",
                timeout: 45000,
                closeOnClick: false,
            });
        });

        this.refreshRuns();
        this.refreshViews();
  }
};
</script>

<style>
    :root {
        --primary: #226597;
        --secondary: #303030;
        --light: #F3F9FB;
        --dark: #113f67;
    }

    .app_left_nav {
        min-width: fit-content;
        border-right: 1px solid #ddd;
    }

    .cutoffText {
        display: inline-block;
        text-overflow: ellipsis;
        overflow: hidden;
        white-space: nowrap;
        flex-grow: 1;
    }
</style>

<!--<style>-->
<!--  @import "https://stackpath.bootstrapcdn.com/bootstrap/4.2.1/css/bootstrap.min.css";-->
<!--  @import "https://cdn.jsdelivr.net/npm/bootstrap-vue@2.0.0-rc.11/dist/bootstrap-vue.css";-->
<!--#app {-->
<!--  width:100%;-->
<!--  background-color: #FFFFFF;-->
<!--}-->
<!--.main{-->
<!--  height: 90vh;-->
<!--  overflow: auto;-->
<!--}-->

<!--.font-text {-->
<!--  font-family: Lato;-->
<!--  font-style: normal;-->
<!--}-->

<!--.font-size-xl {-->
<!--  font-size: 40px;-->
<!--  line-height: 48px;-->
<!--}-->

<!--.font-size-l {-->
<!--  font-size: 24px;-->
<!--  line-height: 29px;-->
<!--}-->

<!--.font-size-m {-->
<!--  font-size: 20px;-->
<!--  line-height: 24px;-->
<!--}-->

<!--.font-size-s {-->
<!--  font-size: 18px;-->
<!--  line-height: 22px;-->
<!--}-->

<!--.font-size-xs {-->
<!--  font-size: 16px;-->
<!--  line-height: 19px;-->
<!--}-->

<!--.font-size-xxs {-->
<!--  font-size: 14px;-->
<!--  line-height: 17px;-->
<!--}-->

<!--.trecx-color-main {-->
<!--  color: #226597;-->
<!--}-->

<!--.trecx-bg-color {-->
<!--  background-color: #2265970D;-->
<!--}-->

<!--.trecx-bg-color-sec {-->
<!--  background-color: #E8E8E8;-->
<!--}-->

<!--.trecx-text-color {-->
<!--  color: #303030;-->
<!--}-->

<!--.trecx-text-sec-color {-->
<!--  color: #989797;-->
<!--}-->

<!--.f-s {-->
<!--  display: flex;-->
<!--  flex-direction: row;-->
<!--  align-items: flex-start;-->
<!--}-->

<!--.pos-r {-->
<!--  position: relative;-->
<!--}-->

<!--.pos-a {-->
<!--  position: absolute;-->
<!--}-->

<!--.align-items-center {-->
<!--  align-items: center;-->
<!--}-->

<!--.align-text-center {-->
<!--  text-align: center;-->
<!--}-->

<!--.cursor-p {-->
<!--  cursor: pointer;-->
<!--}-->

<!--.component-card {-->
<!--  background: #FFFFFF;-->
<!--  border-radius: 10px;-->
<!--  height: 100%;-->
<!--  overflow-y: scroll;-->
<!--  /* direction: ltr; */-->
<!--}-->
<!--.page-view {-->
<!--  background-color: #E5E5E5; -->
<!--  padding: 20px;-->
<!--  width:100%; -->
<!--  height:100%-->
<!--}-->

<!--  /* .nav-pills .nav-link.active, .nav-pills .show>.nav-link {-->
<!--    color: #fff;-->
<!--    background-color: #28a645;-->
<!--  } */-->
<!--  /* .card-title {-->
<!--    font-size: 18px;-->
<!--    font-weight: bold;-->
<!--    color: #5e6b7e;-->
<!--    border-bottom: 2px solid #28a645;-->
<!--    padding-bottom: 5px;-->
<!--  } */-->
<!--  /* .card-text {-->
<!--    font-size: 13px;-->
<!--    color: #5e6b7e;-->
<!--  } */-->
<!--  /* .card-header{-->
<!--    font-size: 14px !important;-->
<!--  }-->
<!--  td{-->
<!--    font-size: 13px !important;-->
<!--  } */-->
<!--.c-nav{-->
<!--  display:none;-->
<!--}-->
<!--.table-header-class {-->
<!--  background-color:rgba(34, 101, 151, 0.1);-->
<!--}-->

<!--.run-thead-class {-->
<!--  background: #E9F0F5;-->
<!--  border-radius: 10px;-->
<!--}-->


<!--/* .table-header-class::before {-->
<!--  content:"";-->
<!--  background-color:rgba(34, 101, 151, 0.1);-->
<!--  position: absolute;-->
<!--  height: 50px;-->
<!--  width: 100%;-->
<!--  left: 0px;-->
<!--} */-->

<!--/* .table-body-tr-class::after {-->
<!--  content: "";-->
<!--  position: absolute;-->
<!--  border: 0.5px solid #D9D9D9;-->
<!--} */-->
<!--</style>-->
