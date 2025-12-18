<template>
    <b-container class="w-100 d-flex flex-column" style="max-width: none; overflow-y: hidden; flex-grow: 1;">
        <b-row class="header row">
            <b-col 
                v-for="(column, j) in columns" v-bind:key="column" 
                class="cell" :class="{ 'justify-content-center': j != 0 }"
            >
                <b-form-checkbox 
                    v-if="j == 0 && sortedItems.some(item => canBeSelected(item)) && canSelectMultiple" 
                    v-model="allAreSelected" class="checkbox" @change="selectAll"
                />
                <span> {{ column[1] }} </span>
                <b-button 
                    variant="link" class="ml-2" v-if="sorters[j] != null"
                    style="float: right; width: 1em; height: 100%; padding-right: 0.15rem" @click="clickSorter(j)"
                    :aria-sort="(sorters[j] != null) ? (sortingMode.index == j) ? sortingMode.mode : 'none' : null"
                ></b-button>
            </b-col>
        </b-row>
        
        <div style="overflow-y: scroll; overflow-x: hidden; flex-grow: 1;">
            <b-row class="run_item row" v-for="item in sortedItems" v-bind:key="item">
                <b-col v-for="(column, j) in columns" v-bind:key="column" class="cell" :class="{ 'justify-content-center': j != 0 }">
                    <b-form-checkbox v-if="j == 0 && canBeSelected(item)" @change="changeIsSelectedAt(item)" v-model="isSelected[item[identifier]]" class="checkbox"/>
                    <slot :name="column[0]" v-bind:item="item"> <span> {{ item[column[0]] }} </span> </slot>
                </b-col>
            </b-row>
        </div>
    </b-container>
</template>

<script>
    export default {
        props: ["items", "columns", "canSelectMultiple", "sorters", "identifier"],
        data() {
            return {
                isSelected: {},
                allAreSelected: false,
                sortingMode: {
                    index: -1,
                    mode: "ascending"
                }, 
                sortedItems: []
            };
        },
        calculated: {

        },
        methods: {
            sortFn(a, b) {
                let tieBreaker = (a.___index - b.___index) / 1000;

                if (this.sortingMode.mode == "descending") 
                    return tieBreaker + this.sorters[this.sortingMode.index](b, a);
                else
                    return tieBreaker + this.sorters[this.sortingMode.index](a, b);
            },
            setSortedItems() {
                let items = this.items.slice();

                // tieBreaker property named ___index to avoid name clashes
                items.forEach((item, i) => item.___index = i);

                if (this.sortingMode.index >= 0)
                    items.sort((a, b) => this.sortFn(a, b));

                this.sortedItems = items;
            },
            updateSelected() {
                this.allAreSelected = Object.values(this.isSelected).every(isSelected => isSelected);
                this.$emit("updateSelected", this.sortedItems.filter(item => this.isSelected[item[this.identifier]]));
            },
            selectAll(allAreSelected) {
                Object.keys(this.isSelected).forEach(key => this.isSelected[key] = allAreSelected);
                this.updateSelected()
            },
            canBeSelected(item) {
                return !("isSelectable" in item) || item.isSelectable;
            },
            changeIsSelectedAt(item) {
                const itemIsSelected = this.isSelected[item[this.identifier]];

                if (!this.canSelectMultiple && itemIsSelected)
                    Object.keys(this.isSelected).forEach(key => this.isSelected[key] = false);

                this.isSelected[item[this.identifier]] = itemIsSelected;

                // console.log({...item}, item[this.identifier]);

                this.updateSelected()
            },
            clickSorter(j) {
                if (this.sorters[j] != null) {
                    this.setSortedItems()

                    if (this.sortingMode.index != j)
                        this.sortingMode = {
                            index: j,
                            mode: "ascending"
                        };
                    else 
                        this.sortingMode.mode = (this.sortingMode.mode == "ascending") ? "descending" : "ascending";

                    this.sortedItems.sort((a,b) => this.sortFn(a, b));
                }
            }
        },
        watch: {
            // Can't have an isSelected watcher beause that disables the @change event (or at least stops changeIsSelectedAt from running)
            isSelected(isSelected) {
                // console.log({...isSelected}, Object.values(isSelected), Object.values(isSelected).every(isSelected => isSelected));
            },
            items() {
                this.setSortedItems()

                this.isSelected = this.sortedItems
                    .filter(item => this.canBeSelected(item))
                    .reduce((acc, item) => {
                        acc[item[this.identifier]] = this.isSelected[item[this.identifier]] || false;

                        return acc;
                    }, {});
                
                this.updateSelected()
            }
        },
        mounted() {
            this.setSortedItems();

            this.isSelected = this.sortedItems
                .filter(item => this.canBeSelected(item))
                .reduce((acc, item) => {
                    acc[item[this.identifier]] = this.isSelected[item[this.identifier]] || false;

                    return acc;
                }, {});
            
            this.updateSelected();

            // console.log("Listview: mounted", this.sortedItems, this.columns, this.isSelected, this.canSelectMultiple)
        }
    }
</script>

<style scoped>
    [aria-sort=none] {
        background-image: url("data:image/svg+xml;charset=utf-8,%3Csvg xmlns='http://www.w3.org/2000/svg' width='101' height='101' preserveAspectRatio='none'%3E%3Cpath opacity='.3' d='M51 1l25 23 24 22H1l25-22zm0 100l25-23 24-22H1l25 22z'/%3E%3C/svg%3E");
    }

    [aria-sort=descending] {
        background-image: url("data:image/svg+xml;charset=utf-8,%3Csvg xmlns='http://www.w3.org/2000/svg' width='101' height='101' preserveAspectRatio='none'%3E%3Cpath opacity='.3' d='M51 1l25 23 24 22H1l25-22z'/%3E%3Cpath d='M51 101l25-23 24-22H1l25 22z'/%3E%3C/svg%3E");
    }
    
    [aria-sort=ascending] {
        background-image: url("data:image/svg+xml;charset=utf-8,%3Csvg xmlns='http://www.w3.org/2000/svg' width='101' height='101' preserveAspectRatio='none'%3E%3Cpath d='M51 1l25 23 24 22H1l25-22z'/%3E%3Cpath opacity='.3' d='M51 101l25-23 24-22H1l25 22z'/%3E%3C/svg%3E");
    }

    [aria-sort] {
        background-position: right 0.15rem center;
        background-repeat: no-repeat;
        background-size: 0.65em 1em;
    }

    .header {
        background: var(--light);
        color: var(--dark);
        font-weight: bold;
        line-height: 22px;
    }

    .row {
        padding: 15px 20px;
        font-size: 16px;
    }

    .run_item {
        line-height: 20px;
        border-bottom: 1px solid #ddd;
    }

    .cell {
        display: flex;
        flex-direction: row;
        align-items: center;
    }

    .checkbox {
        padding-right: 20px;
    }
</style>