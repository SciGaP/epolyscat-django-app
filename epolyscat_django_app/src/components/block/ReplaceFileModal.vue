<template>
    <b-modal title="Please Confirm" ref="modal" @hide="hidden()" centered>
        {{ filename }} has already been uploaded, do you want to replace it?
        <template #modal-footer="{ ok, cancel }">
            <b-button variant="danger" @click="ok(); resolve('yes')">Yes</b-button>
            <b-button variant="outline-danger" @click="ok(); resolve('all')">Replace all</b-button>
            <b-button variant="dark" @click="cancel(); resolve('no')">No</b-button>
            <b-button variant="outline-dark" @click="cancel(); resolve('none')">Replace none</b-button>
        </template>
    </b-modal>
</template>

<script>
export default {
    props: [],
    data() {
        return {
            filename: null,
            resolve: null
        }
    },
    methods: {
        hidden() {
            // this.resolve('none');
        },
        async ask(filename) {
            this.$refs.modal.show();
            this.filename = filename

            return await (new Promise(resolve => {
                this.resolve = resolve;
            }))
        }
    }
}
</script>