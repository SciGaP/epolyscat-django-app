<template>
  <div class="w-100">

    <div class="w-100 left-nav-menu">
      <div class="w-100">
        <router-link :to="`/`" v-slot="{ href, route, navigate, isActive,isExactActive }">
          <b-link :class="{active: isExactActive}" :href="href" @click="navigate">
            <svg class="d-inline-block" width="22" height="23" viewBox="0 0 22 23" fill="none"
                 xmlns="http://www.w3.org/2000/svg">
              <path
                  d="M20.7862 10.1765L11.487 1.19805C11.216 0.933983 10.784 0.933983 10.5131 1.19805L1.2142 10.1761C1.0773 10.3083 1 10.4903 1 10.6806C1 11.0679 1.31394 11.3818 1.70114 11.3818H2.9369V21.1977C2.9369 21.585 3.25079 21.8989 3.63804 21.8989H8.54602C8.93327 21.8989 9.24716 21.585 9.24716 21.1977V15.5886H12.7529V21.1977C12.7529 21.585 13.0668 21.8989 13.454 21.8989H18.362C18.7492 21.8989 19.0631 21.585 19.0631 21.1977V11.3818H20.2989C20.489 11.3818 20.6709 11.3046 20.803 11.1679C21.0721 10.8895 21.0646 10.4456 20.7862 10.1765ZM19.0631 10.6806C18.6759 10.6806 18.362 10.9945 18.362 11.3818V21.1977H13.454V15.5886C13.454 15.2014 13.1401 14.8875 12.7529 14.8875H9.24716C8.85991 14.8875 8.54602 15.2014 8.54602 15.5886V21.1977H3.63804V11.3818C3.63804 10.9945 3.32415 10.6806 2.9369 10.6806H1.70114L11 1.70252L20.3024 10.6806H19.0631Z"
                  fill="#226597" stroke="#226597"/>
            </svg>
            Home
          </b-link>
        </router-link>
      </div>
    </div>

    <div class="w-100 left-nav-menu">
      <div class="w-100">
        <router-link :to="tutorialsLink" v-slot="{ href, route, navigate, isActive,isExactActive }">
          <b-overlay :show="!tutorialsLink">
            <b-link :class="{active: isExactActive}" :href="href" @click="navigate">
              <svg class="d-inline-block" width="20" height="20" viewBox="0 0 20 20" fill="none"
                   xmlns="http://www.w3.org/2000/svg">
                <path d="M20 20H0V4.80005H20V20ZM0.8 19.2H19.2V5.60005H0.8V19.2Z" fill="black"/>
                <path d="M16.3996 1.6001H3.59961V2.4001H16.3996V1.6001Z" fill="black"/>
                <path d="M14.7997 0H5.19971V0.8H14.7997V0Z" fill="black"/>
                <path d="M17.9995 3.2002H1.99951V4.0002H17.9995V3.2002Z" fill="black"/>
              </svg>
              <div class="d-inline-block" style="flex: 1;padding-left:10px">Tutorials</div>
            </b-link>
            <b-icon icon="caret-right-fill" scale="0.5"></b-icon>
          </b-overlay>
        </router-link>

      </div>
    </div>

    <div class="w-100 left-nav-menu">
      <div class="w-100" v-b-toggle.collapse-recent-experiments-list>
        <svg class="d-inline-block" width="22" height="20" viewBox="0 0 22 20" fill="none"
             xmlns="http://www.w3.org/2000/svg">
          <path
              d="M21 9.95731C21 8.49887 19.4953 7.22559 17.1236 6.42C17.1829 6.12286 17.2317 5.83079 17.267 5.54629C17.5294 3.42934 17.0795 1.92016 16 1.297C14.9204 0.673637 13.3884 1.03875 11.6865 2.32438C11.4573 2.49754 11.2282 2.68637 11 2.88684C10.7718 2.68637 10.5427 2.49754 10.3135 2.32438C8.61176 1.03875 7.07965 0.67434 6 1.297C4.92055 1.9202 4.47059 3.42934 4.73301 5.54629C4.76832 5.83075 4.81715 6.12282 4.87645 6.42C2.50473 7.22559 1 8.49887 1 9.95731C1 11.4155 2.50422 12.6886 4.87523 13.4942C4.38789 15.9495 4.73742 17.8886 6 18.6176C6.34414 18.8161 6.73441 18.9145 7.16141 18.9145C8.07359 18.9145 9.15414 18.4659 10.3135 17.5902C10.5427 17.417 10.7718 17.2282 11 17.0277C11.2282 17.2282 11.4573 17.417 11.6865 17.5902C12.8459 18.466 13.9263 18.9145 14.8386 18.9145C15.2656 18.9145 15.6558 18.8162 16 18.6176C17.2626 17.8886 17.6121 15.9495 17.1248 13.4941C19.4958 12.6886 21 11.4155 21 9.95731ZM12.133 2.91555C13.1532 2.14497 14.0945 1.74313 14.8384 1.74313C15.1363 1.74313 15.4027 1.8075 15.6296 1.93864C16.4234 2.39672 16.7521 3.67856 16.5319 5.45516C16.5018 5.69805 16.462 5.94661 16.4132 6.19981C15.5271 5.95059 14.5425 5.76207 13.4865 5.64614C12.8682 4.80473 12.2142 4.0445 11.546 3.39207C11.7413 3.22274 11.9372 3.06348 12.133 2.91555ZM6.44883 10.8682C6.66566 11.3065 6.8998 11.7446 7.1509 12.1795C7.40219 12.6146 7.66461 13.0367 7.9359 13.4438C7.15773 13.3306 6.42434 13.1761 5.75309 12.9872C5.92512 12.3113 6.1577 11.5991 6.44883 10.8682ZM6.44652 9.0502C6.15703 8.32637 5.92223 7.61387 5.74945 6.92844C6.42172 6.73911 7.15629 6.58426 7.93594 6.47086C7.66465 6.87789 7.40223 7.3 7.15094 7.73512C6.89695 8.17493 6.6625 8.61411 6.44652 9.0502ZM7.79254 11.8091C7.43809 11.1953 7.11953 10.5759 6.8373 9.96008C7.11758 9.34536 7.43555 8.72375 7.79254 8.10543C8.14723 7.49114 8.52453 6.90516 8.91707 6.3525C9.58906 6.28793 10.2866 6.2536 11 6.2536C11.7134 6.2536 12.4109 6.28797 13.0829 6.3525C13.4755 6.90516 13.8528 7.49114 14.2075 8.10543C14.5645 8.72375 14.8825 9.34536 15.1627 9.96008C14.8805 10.5759 14.5619 11.1953 14.2075 11.8091C13.8528 12.4234 13.4755 13.0094 13.0829 13.5621C12.4109 13.6266 11.7134 13.661 11 13.661C10.2866 13.661 9.58906 13.6266 8.91707 13.5621C8.52453 13.0094 8.14723 12.4234 7.79254 11.8091ZM14.8491 12.1795C15.1002 11.7446 15.3344 11.3065 15.5512 10.8683C15.8423 11.5991 16.0748 12.3113 16.2469 12.9872C15.5757 13.1761 14.8423 13.3306 14.0641 13.4438C14.3354 13.0367 14.5979 12.6147 14.8491 12.1795ZM14.8491 7.73508C14.5979 7.29996 14.3354 6.87789 14.0641 6.47082C14.8437 6.58422 15.5784 6.73907 16.2506 6.9284C16.0779 7.61383 15.843 8.32633 15.5535 9.05016C15.3376 8.61411 15.103 8.17493 14.8491 7.73508ZM11 3.89368C11.5079 4.38641 12.0082 4.94665 12.491 5.56008C12.0055 5.52911 11.5078 5.51286 11 5.51286C10.4922 5.51286 9.99461 5.52911 9.50906 5.56004C9.9918 4.94664 10.4921 4.38641 11 3.89368ZM6.37039 1.93864C6.59734 1.80754 6.86355 1.74313 7.1616 1.74313C7.90559 1.74313 8.84707 2.14497 9.86703 2.91555C10.0629 3.06348 10.2587 3.2227 10.4541 3.39207C9.7859 4.04446 9.13195 4.80469 8.51359 5.64614C7.45766 5.76207 6.47293 5.95059 5.58695 6.19981C5.53816 5.94661 5.49832 5.69805 5.4682 5.45516C5.24785 3.67856 5.57664 2.39672 6.37039 1.93864ZM1.74074 9.95731C1.74074 8.8479 3.03605 7.83383 5.04164 7.14829C5.2725 8.05301 5.60383 8.99926 6.02328 9.95532C5.59566 10.9279 5.26797 11.8741 5.04082 12.7661C3.0357 12.0806 1.74074 11.0666 1.74074 9.95731ZM9.86699 16.999C8.43848 18.0783 7.16461 18.4346 6.37035 17.9759C5.41008 17.4216 5.18016 15.7934 5.58848 13.7152C6.4741 13.9642 7.45824 14.1526 8.51355 14.2684C9.13191 15.1099 9.78582 15.8701 10.4541 16.5225C10.2587 16.6919 10.0629 16.8511 9.86699 16.999ZM11 16.0209C10.4921 15.5282 9.9918 14.968 9.50906 14.3545C9.99461 14.3855 10.4922 14.4017 11 14.4017C11.5078 14.4017 12.0055 14.3855 12.491 14.3545C12.0082 14.968 11.5079 15.5282 11 16.0209ZM15.6296 17.9759C14.8364 18.4349 13.5615 18.0785 12.133 16.999C11.9371 16.8511 11.7413 16.6919 11.5459 16.5225C12.2142 15.8701 12.8681 15.1099 13.4864 14.2684C14.5418 14.1526 15.5259 13.9642 16.4115 13.7152C16.8198 15.7934 16.59 17.4216 15.6296 17.9759ZM16.9592 12.7661C16.732 11.8741 16.4043 10.9279 15.9768 9.95532C16.3962 8.99926 16.7275 8.05301 16.9584 7.14829C18.964 7.83383 20.2593 8.8479 20.2593 9.95731C20.2593 11.0666 18.9643 12.0806 16.9592 12.7661Z"
              fill="black" stroke="black" stroke-width="0.3"/>
        </svg>
        <div class="d-inline-block" style="flex: 1;padding-left:10px">Recent Experiments</div>
        <div class="d-inline-block">
          <span class="when-open"><b-icon icon="caret-down-fill" scale="0.5"></b-icon></span>
          <span class="when-closed"><b-icon icon="caret-right-fill" scale="0.5"></b-icon></span>
        </div>
      </div>
      <table-overlay-info rows="5" columns="1" :data="experiments" empty-label="No experiments">
        <b-collapse visible id="collapse-recent-experiments-list" accordion="left-nav-sub-list" class="w-100 pl-4 pt-2">
          <ol class="w-100 overflow-auto" style="max-height: 400px;">
            <li v-for="(experiment) in experiments" :key="experiment.experimentId">
              <router-link :to="`/runs/?experimentId=${experiment.experimentId}`"
                           v-slot="{ href, route, navigate, isActive,isExactActive }">
                <b-link :class="{active: isExactActive}" :href="href" @click="navigate">
                  <div style="display: flex; flex-direction: row;">
                    <div class="overflow-auto" style="flex: 1;">{{ experiment.name }}</div>
                    <b-badge>{{ experiment.activeRunCount }} runs</b-badge>
                  </div>
                </b-link>
              </router-link>
            </li>
          </ol>
          <router-link :to="`/experiments`"
                       v-slot="{ href, route, navigate, isActive,isExactActive }">
            <b-link :href="href" @click="navigate">View all
              <b-icon icon="chevron-right" scale="0.5"></b-icon>
            </b-link>
          </router-link>
        </b-collapse>
      </table-overlay-info>
    </div>

    <div class="w-100 left-nav-menu">
      <div class="w-100" v-b-toggle.collapse-views-list>
        <svg class="d-inline-block" width="20" height="20" viewBox="0 0 20 20" fill="none"
             xmlns="http://www.w3.org/2000/svg">
          <path
              d="M18.9102 10.3061L19.0977 10.0001L18.9102 9.69381C18.8787 9.64269 18.1243 8.42503 16.6587 7.20477C14.7094 5.58198 12.413 4.72412 10.0175 4.72412C7.62262 4.72412 5.32053 5.58152 3.36008 7.20355C1.88608 8.42319 1.12314 9.63994 1.0914 9.69121L0.899902 10.0001L1.0914 10.3089C1.12314 10.3602 1.88608 11.5769 3.36008 12.7966C5.32053 14.4186 7.62262 15.276 10.0175 15.276C12.413 15.276 14.7094 14.4181 16.6587 12.7953C18.1243 11.5751 18.8787 10.3574 18.9102 10.3061ZM15.8778 11.9204C14.1242 13.3694 12.1525 14.1039 10.0175 14.1039C7.88171 14.1039 5.90356 13.3689 4.13812 11.919C3.22504 11.1692 2.60599 10.4099 2.3028 10.0001C2.60614 9.58989 3.22519 8.83076 4.13812 8.08094C5.90356 6.63118 7.88156 5.89601 10.0175 5.89601C12.1525 5.89601 14.1241 6.63073 15.8778 8.07956C16.7867 8.83061 17.4005 9.59096 17.6995 10.0001C17.4001 10.4093 16.7866 11.1697 15.8778 11.9204Z"
              fill="black"/>
          <path
              d="M9.99982 6.479C8.05844 6.479 6.479 8.05846 6.479 9.99985C6.479 11.9412 8.05844 13.5207 9.99982 13.5207C11.9412 13.5207 13.5206 11.9412 13.5206 9.99985C13.5206 8.05846 11.9412 6.479 9.99982 6.479ZM9.99982 12.3488C8.70465 12.3488 7.65088 11.295 7.65088 9.99985C7.65088 8.70467 8.70465 7.65089 9.99982 7.65089C11.295 7.65089 12.3488 8.70467 12.3488 9.99985C12.3488 11.295 11.295 12.3488 9.99982 12.3488Z"
              fill="black"/>
          <path d="M1.17188 1.17189H3.61816V0H0V3.57426H1.17188V1.17189Z" fill="black"/>
          <path d="M16.3818 0V1.17189H18.8281V3.57426H20V0H16.3818Z" fill="black"/>
          <path d="M1.17188 16.4255H0V19.9998H3.61816V18.8279H1.17188V16.4255Z" fill="black"/>
          <path d="M18.8281 18.8279H16.3818V19.9998H20V16.4255H18.8281V18.8279Z" fill="black"/>
        </svg>
        <div class="d-inline-block" style="flex: 1;padding-left:10px">Views</div>
        <div class="d-inline-block">
          <span class="when-open"><b-icon icon="caret-down-fill" scale="0.5"></b-icon></span>
          <span class="when-closed"><b-icon icon="caret-right-fill" scale="0.5"></b-icon></span>
        </div>
      </div>
      <b-collapse visible id="collapse-views-list" accordion="left-nav-sub-list" class="w-100 pl-4 pt-2 overflow-auto" style="max-height: 400px;">
        <ol class="w-100">
          <li v-for="(view) in views" :key="view.viewId">
            <router-link :to="`/runs/?viewId=${view.viewId}`"
                         v-slot="{ href, route, navigate, isActive,isExactActive }">
              <b-link :class="{active: isExactActive}" :href="href" @click="navigate">
                <div style="display: flex; flex-direction: row;">
                  <div class="overflow-auto" style="flex: 1;">{{ view.name }}</div>
                  <b-badge>{{ view.activeRunCount }} runs</b-badge>
                </div>
              </b-link>
            </router-link>
          </li>
        </ol>
        <router-link :to="`/views`"
                     v-slot="{ href, route, navigate, isActive,isExactActive }">
          <b-link :href="href" @click="navigate">View all
            <b-icon icon="chevron-right" scale="0.5"></b-icon>
          </b-link>
        </router-link>
      </b-collapse>
    </div>

    <div class="w-100 left-nav-menu">
      <div class="w-100">

        <!--        <div class="d-inline-block" style="padding-left:10px"><span>All runs</span></div>-->
        <router-link :to="`/runs`" v-slot="{ href, route, navigate, isActive,isExactActive }">
          <b-link :class="{active: isExactActive}" :href="href" @click="navigate">
            <svg class="d-inline-block" width="20" height="20" viewBox="0 0 20 20" fill="none"
                 xmlns="http://www.w3.org/2000/svg">
              <path d="M20 20H0V4.80005H20V20ZM0.8 19.2H19.2V5.60005H0.8V19.2Z" fill="black"/>
              <path d="M16.3996 1.6001H3.59961V2.4001H16.3996V1.6001Z" fill="black"/>
              <path d="M14.7997 0H5.19971V0.8H14.7997V0Z" fill="black"/>
              <path d="M17.9995 3.2002H1.99951V4.0002H17.9995V3.2002Z" fill="black"/>
            </svg>
            <div class="d-inline-block" style="flex: 1;padding-left:10px">All runs</div>
          </b-link>
        </router-link>
        <b-icon icon="caret-right-fill" scale="0.5"></b-icon>
      </div>
    </div>

  </div>
</template>


<script>
import store from "@/store";
import TableOverlayInfo from "@/components/overlay/table-overlay-info";

export default {
  name: 'AppLeftNav',
  components: {TableOverlayInfo},
  store: store,
  computed: {
    experiments() {
      return this.$store.getters["experiment/getExperiments"]();
    },
    views() {
      return this.$store.getters["view/getViews"]();
    },
    tutorialViews() {
      return this.$store.getters["view/getViews"]({tutorials: true});
    },
    tutorialsLink() {
      if (this.tutorialViews && this.tutorialViews.length > 0) {
        return `/runs/?viewId=${this.tutorialViews[0].viewId}`;
      } else {
        return null;
      }
    }
  },
  async mounted() {
    await this.$store.dispatch("experiment/fetchExperiments");
    await this.$store.dispatch("view/fetchViews");
    await this.$store.dispatch("view/fetchViews", {tutorials: true});
  }
}
</script>

<style scoped>
.left-nav-menu {
  padding: 10px;
}

.left-nav-menu ol {
  padding: 0px;
  margin: 0px;
  list-style: none;
  border-left: 2px solid rgba(35, 122, 187, .4);
  display: inline-block;
}

.left-nav-menu ol li {
  padding: 10px;
  line-height: 17px;
}

/*.app-left-nav {*/
/*  width: 276px !important;*/
/*  min-width: 276px !important;*/
/*  height: 100%;*/
/*  box-sizing: border-box;*/
/*  border: 1px solid #E5E5E5;*/
/*}*/

/*.nav-content {*/
/*  margin-top: 59px;*/
/*  margin-left: 32px;*/
/*  margin-right: 26px;*/
/*  padding: 0px;*/
/*}*/

/*li {*/
/*  list-style-type: none;*/
/*  margin-bottom: 22px !important;*/
/*}*/

/*.run-badges {*/
/*  float: right;*/
/*}*/

/*.badge {*/
/*  width: 55px !important;*/
/*  height: 30px !important;*/
/*  font-size: 13.3333px;*/
/*  line-height: 16px;*/
/*  text-align: center;*/
/*  float: right;*/
/*  border-radius: 5px;*/
/*}*/

/*.sub-list ol {*/
/*  position: relative;*/
/*  padding-left: 40px !important;*/
/*}*/

/*.sub-list ol li {*/
/*  position: relative;*/
/*  list-style: none;*/
/*  padding-top: 14px !important;*/
/*  margin-bottom: 0px !important;*/
/*}*/

/*.sub-list ol:before {*/
/*  content: " ";*/
/*  display: block;*/
/*  position: absolute;*/
/*  top: 14px !important;*/
/*  bottom: 0;*/
/*  left: 32px !important;*/
/*  border-left: 2px solid rgba(35, 122, 187, 0.4);*/
/*  border-radius: 2px;*/
/*}*/


</style>