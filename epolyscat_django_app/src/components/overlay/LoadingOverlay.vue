<template>
    <b-overlay :show="messages.length != 0" style="height: fit-content" class="w-100 full_width">
        <slot></slot>
        <template #overlay>
            <div class="d-flex flex-row align-items-center justify-content-center">
                <b-spinner label="Spinning" style="flex-shrink: 0"/>
                <h5 class="cutoffText" style="margin-left: 10px; margin-bottom: 0; flex-grow: 0">
                    {{ messages.join(", ") }}
                </h5>
            </div>
        </template>
    </b-overlay>
</template>

<script>
import store from '@/store';

export default {
    props: ["name"],
    store,
    data() {
        return {
            loadingEvents: []
        };
    },
    computed: {
        messages() {
            return this.$store.getters["loading/getMessages"](this.name);
        }
    }
}
</script>

<style>
.full_width > .b-overlay > * {
    width: 100%!important;
    padding: 20px;
}
</style>