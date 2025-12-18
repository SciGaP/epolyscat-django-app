<template>
    <div class="d-flex flex-column container w-100">
        <h5 style="word-break: break-word;">{{ name }}</h5>
        <hr style="width: 100%; background: white">
        <span v-if="error.message"><b>Error message:</b> {{ error.message }}</span>
        <span v-if="errorExplaination != null">
            <b>Explanation:</b> {{ errorExplaination }}
        </span>
        <span v-if="needsAuthentication">
            You might need to login again, you can do that 
            <b-link variant="link" @click="$router.go()" class="text-light" style="text-decoration: underline">here</b-link>
            or you can <b-link variant="link" @click="loginSeperateTab()" class="text-light" style="text-decoration: underline">login in a seperate tab</b-link>.
        </span>
        <b-button 
            variant="link" size="sm" @click="showMore = !showMore"
            :aria-expanded="showMore ? 'true' : 'false'"
            aria-controls="moreInfo" class="text-light mb-2 mt-4 p-0" style="width: fit-content;"
        >
            {{ (showMore) ? "Show Less" : "Show More" }}
        </b-button>
        <b-collapse 
            v-model="showMore" id="moreInfo" class="m-0 p-0" 
            style="white-space: pre; overflow: scroll; max-height: 600px; width: 450px;"
        > {{ error.stack }} </b-collapse>
    </div>
</template>

<script>
import router from '@/router';

export default {
    props: [ "name", "error" ],
    router,
    data() {
        return {
            showMore: false
        }
    },
    computed: {
        needsAuthentication() {
            const err = this.error.message;

            return this.includesErrorCode(err, "403") ||
                this.includesErrorCode(err, "401");
        },
        errorExplaination() {
            const err = this.error.message;

            if (this.includesErrorCode(err, "500"))
                return "Server error";

            return null;
        }
    },
    methods: {
        includesErrorCode(errorMessage, code) {
            return errorMessage.split("id:")[0].split("name:")[0].includes(code);
        },
        loginSeperateTab() {
            open(window.location.origin + "/epolyscat_django_app/home");
        }
    },
    mounted() {
        console.log(this.name, this.error);
    }
}
</script>
