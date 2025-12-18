<template>
    <div class="tableItem">
        <span v-if="'name' in item && includeName" style="white-space: nowrap; margin-right: 10px"> 
            {{ item.name }} 
        </span>
        <b-form-select 
            v-if="item.type == 'select'" :value="item.data"
            @input="set($event)" :options="item.options"
            :disabled="viewing" size="sm"
            :title="('tooltip' in item) ? item.tooltip : null"
            v-b-tooltip.hover.right
        />
        <b-form-group class="w-100 m-0" v-else>
            <b-form-input 
                size="sm"
                :value="item.data"
                @input="set($event)"
                :aria-describedby="`input-help-${id} input-feedback-${id}`"
                :type="('type' in item) ? item.type : 'text'"
                :placeholder="('placeholder' in item) ? item.placeholder : ''"
                :min="('min' in item) ? item.min : -Infinity"
                :max="('max' in item) ? item.max : Infinity"
                :state="('issues' in item && !viewing) ? issues.length == 0 : null"
                :disabled="viewing"
                :step="('step' in item) ? item.step : 1"
                :title="('tooltip' in item) ? item.tooltip : null"
                v-b-tooltip.hover.right
            />
            <b-form-text v-if="'description' in item && item.description">{{item.description}}</b-form-text>
            <template v-slot:invalid-feedback>
                <div 
                    v-for="issue in issues" 
                    :key="issue" :id="`input-feedback-${id}`"
                >{{ issue }}</div>
            </template>
            <template v-slot:description v-if="'help' in item">
                {{ item.help }}
            </template>
        </b-form-group>
    </div>
</template>

<script>
export default {
    props: ["item", "includeName", "viewing"],
    data () {
        return {
            id: null
        }
    }, 
    computed: {
        issues() {
            return ("issues" in this.item) ? this.item.issues(this.item.data) : [];
        }
    },
    methods: {
        set(data) {
            this.$emit("changed")
            this.item.data = data;
            
            if (this.issues.length == 0)
                this.item.set(data);
        }
    },
    mounted () {
        // console.log("Data: ", this.item.data, JSON.stringify(this.item), this.item);
        this.id = this._uid;
        // this.item.data = ("data" in this.item && this.item.data != null) ? this.item.data: "";
    }
}
</script>

<style scoped>
.tableItem {
    display: flex;
    flex-direction: column;
    width: 100%
}
</style>