import Vue from 'vue';
import { BootstrapVue, BootstrapVueIcons } from 'bootstrap-vue';
// import userPlugin from "./components/common/UserPlugin.js";
// import StarRating from 'vue-star-rating';

// import { FontAwesomeIcon } from '@fortawesome/vue-fontawesome';
import router from './router';
// import store from "./store";
import App from './App.vue';
import VueFilterDateFormat from 'vue-filter-date-format';
import Toast from "vue-toastification";
import "vue-toastification/dist/index.css";

// Vue.component('font-awesome-icon', FontAwesomeIcon);
// Vue.component('star-rating', StarRating);

Vue.config.productionTip = false;
Vue.use(BootstrapVue);
Vue.use(BootstrapVueIcons);
// Vue.use(userPlugin);
Vue.use(VueFilterDateFormat);

//toast options
// const toasterOptions = {
//   // You can set your default options here
//   timeout: 5000,
// };
Vue.use(Toast, {
    maxToasts: 5,
    filterBeforeCreate: (newToast, toasts) => {
        return (toasts.length > 0 && toasts.some(toast =>
            toast.content.props.name == newToast.content.props.name &&
            toast.content.props.error.message == newToast.content.props.error.message
        )) ? false : newToast;
    }
});


new Vue({
  // store,
  router,
  render: h => h(App),
}).$mount('#app');
