<template>
    <div class="d-flex flex-column align-items-center">
        <span class="mb-1">{{ label }}</span>
        <div class="group_box d-flex flex-row">
            <div 
                class="group_box_part btn-group-toggle" data-toggle="buttons" 
                v-for="(buttons, i) in buttonsList" v-bind:key="buttons" 
                :style=" { 'grid-template-rows': `repeat(${buttons.perColumn}, 1fr)` } "
            >
                <label v-for="button in buttons" v-bind:key="button"
                    v-b-tooltip.noninteractive.hover.top="(button.checked == undefined || button.checked) ? false : 'Still needs to be uploaded'"
                    class="btn icon-spacing" :class="{
                        'btn-primary': path[i + offset] == button.key, 
                        'btn-light': path[i + offset] != button.key
                    }" @click="$emit('setPath', {index: i + offset, key: button.key})">
                    <input type="radio" :name="'keys' + i" :id="button.key" autocomplete="off"> 
                    <span style="margin: 0 auto">{{ button.key }}</span>
                    <template v-if="button.checked != undefined">
                        <b-icon icon="file-earmark-check" v-if="button.checked" />
                        <b-icon icon="file-earmark-plus"     v-else />
                    </template>
                    <a v-if="button.delete" @click="button.delete()">
                        <b-icon icon="trash" size="sm" />
                    </a>
                </label>
            </div>
        </div>
    </div>
</template>

<script>


export default {
    props: ["buttonsList", "offset", "label"],
    computed: {
        path() {
            return this.$store.getters["input/getPath"];
        }
    },
}
</script>

<style scoped>
.group_box {
    border: 1px solid #ddd;
    border-radius: 8px;
    overflow: hidden;
    width: fit-content;
}

.group_box_part {
    overflow-y: scroll;
    outline: 0.5px solid #ddd;
    /* margin-left: 1px; */
    display: grid;
    grid-auto-flow: column;
    padding: 10px;
    max-height: 240px;
}

.group_box_part > * {
    margin: 5px 5px;
}

.icon-spacing {
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    align-items: center;
}

.icon-spacing > * {
    margin-left: 15px;
}
</style>