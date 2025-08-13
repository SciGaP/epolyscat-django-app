<template>
  <div class="w-100" v-if="inputTable">
    <b-tabs pills content-class="w-100 overflow-auto" nav-class="pages-nav-pills">
      <b-tab :title="page.name" v-for="(page, pageIndex) in inputTable.pages" :key="pageIndex">
        <hr class="w-100"/>
        <b-tabs pills content-class="w-100 overflow-auto" nav-class="sections-nav-pills">
          <div class="w-100">
            <b-tab :title="`${section.category}`" v-for="(section, sectionIndex) in page.sections"
                   :key="`${pageIndex}-${sectionIndex}`">
              <div class="overflow-auto bg-light p-2">
                <b-table-simple class="w-100 overflow-auto">
                  <b-thead>
                    <b-tr>
                      <b-th v-for="(header, headerIndex) in section.lines[0].items"
                            :key="`${pageIndex}-${sectionIndex}-th-${headerIndex}`">
                        {{ header.name }}
                      </b-th>
                      <b-th></b-th>
                    </b-tr>
                  </b-thead>
                  <b-tbody>
                    <b-tr v-for="(line, lineIndex) in section.lines"
                          :key="`${pageIndex}-${sectionIndex}-tr-${lineIndex}`">
                      <b-td v-for="(object, objectIndex) in line.items" style="min-width: 120px;"
                            :key="`${pageIndex}-${sectionIndex}-tr-${lineIndex}-td-${objectIndex}`">
                        <b-form-input
                            :id="`${pageIndex}-${sectionIndex}-tr-${lineIndex}-td-${objectIndex}-input`"
                            size="sm" :readonly="readOnly"
                            v-model="object.value" v-on:change="onChange"
                            :title="object.documentation"
                            :placeholder="object.defaultValue"/>
                      </b-td>
                      <b-td>
                        <div class="d-flex flex-row">
                          <b-button size="sm" variant="link" v-on:click="deleteLine(section, lineIndex)"
                                    v-b-tooltip.hover title="Delete" :disabled="readOnly || section.lines.length < 2">
                            <b-icon icon="trash"/>
                          </b-button>
                        </div>
                      </b-td>
                    </b-tr>
                  </b-tbody>
                </b-table-simple>
                <b-button variant="link" v-on:click="addNewLine(section)" :disabled="readOnly">Add new row</b-button>
              </div>
            </b-tab>
          </div>
        </b-tabs>
      </b-tab>
    </b-tabs>
  </div>
</template>

<script>
export default {
  name: "RunInputTable",
  props: {
    "inputTable": {},
    "readOnly": {
      default: false
    }
  },
  data() {
    return {}
  },
  methods: {
    addNewLine(section) {
      if (section.lines.length > 0) {
        section.lines.push({
          items: section.lines[section.lines.length - 1].items.map(item => {
            return {...item, value: ""}
          })
        })
      }

      this.onChange();
    },
    deleteLine(section, lineIndex) {
      // Delete avoid if there's only one line.
      if (section.lines.length > 1) {
        section.lines = section.lines.filter((line, i) => i !== lineIndex);
      }

      this.onChange();
    },
    onChange() {
      this.$emit("change", this.inputTable);
    }
  }
}
</script>

<style lang="scss">
@import "../../../styles";

ul.pages-nav-pills {
  //margin-left: 150px;
}

ul.pages-nav-pills,
ul.sections-nav-pills {
}

ul.pages-nav-pills li.nav-item,
ul.sections-nav-pills li.nav-item {
}

ul.pages-nav-pills li.nav-item a.nav-link,
ul.sections-nav-pills li.nav-item a.nav-link {
  border-radius: 20px;
  padding: 5px 30px;
  margin-right: 10px;
  margin-bottom: 10px;
  background-color: $light;
}

ul.pages-nav-pills li.nav-item a.nav-link.active {
  background-color: $primary;
  color: $light;
}

ul.pages-nav-pills li.nav-item a.nav-link {
  background-color: $light;
  color: $primary;
}

ul.sections-nav-pills li.nav-item a.nav-link.active {
  background-color: $info;
  color: $light;
}

ul.sections-nav-pills li.nav-item a.nav-link {
  background-color: $light;
  color: $info;
}
</style>
