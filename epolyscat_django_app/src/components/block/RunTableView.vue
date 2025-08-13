<template>
  <div class="run-table-view">
    <b-table
      :items="paramsList"
      :fields="fields"
      borderless=false
      thead-class="run-thead-class"
    >
      <template v-slot:cell(name)="row">
        <b-form-group>
          <input type="text" v-model="row.item.name" />
        </b-form-group>
      </template>
      <template v-slot:cell(lower_end)="row">
        <b-form-group>
          <input type="text" v-model="row.item.lower_end" />
        </b-form-group>
      </template>
      <template v-slot:cell(upper_end)="row">
        <b-form-group>
          <input type="text" v-model="row.item.upper_end" />
        </b-form-group>
      </template>
      <template v-slot:cell(order)="row">
        <b-form-group>
          <input type="text" v-model="row.item.order" />
        </b-form-group>
      </template>
      <template v-slot:cell(ncofficients)="row">
        <b-form-group>
          <input type="text" v-model="row.item.ncofficients" />
        </b-form-group>
      </template>
      <template v-slot:cell(axis)="row">
        <span class="font-size-s trecx-text-color">{{ row.item.axis }}</span>
      </template>
      <template #head(axis)="row">
        <div style="min-width: 30px !important;"><span class="font-text font-size-m trecx-color-main"><b>{{ row.field.label }}</b></span></div>
      </template>
      <template #head(lower_end)="row">
        <div><span class="font-text font-size-m trecx-color-main"><b>{{ row.field.label }}</b></span></div>
      </template>
      <template #head(upper_end)="row">
        <div><span class="font-text font-size-m trecx-color-main"><b>{{ row.field.label }}</b></span></div>
      </template>
      <template #head(order)="row">
        <div><span class="font-text font-size-m trecx-color-main"><b>{{ row.field.label }}</b></span></div>
      </template>
      <template #head(name)="row">
        <div><span class="font-text font-size-m trecx-color-main"><b>{{ row.field.label }}</b></span></div>
      </template>
      <template #head(ncofficients)="row">
        <div><span class="font-text font-size-m trecx-color-main"><b>{{ row.field.label }}</b></span></div>
      </template>
    </b-table>

  </div>
</template>

<script>
export default {
  name: 'RunTableView',
  data() {
    return {
      fields: [
        {
          key: 'axis',
          label: "Axis",
        },
        {
          key: 'name',
          label: 'Name',
        },
        {
          key: 'ncofficients',
          label: 'N-coefficients',
        },
        {
          key: 'lower_end',
          label: 'Lower end'
        },
        {
          key: 'upper_end',
          label: 'Upper end'
        },
        {
          key: 'order',
          label: 'Order'
        }
      ],
      paramsList: [
        { axis:'[hide]', name: '', ncofficients: '', lower_end: '', upper_end: '', order: ''},
        { axis:'[hide]', name: '', ncofficients: '', lower_end: '', upper_end: '', order: ''},
        { axis:'[hide]', name: '', ncofficients: '', lower_end: '', upper_end: '', order: ''},
        { axis:'[hide]', name: '', ncofficients: '', lower_end: '', upper_end: '', order: ''},
        { axis:'[hide]', name: '', ncofficients: '', lower_end: '', upper_end: '', order: ''},
        { axis:'[hide]', name: '', ncofficients: '', lower_end: '', upper_end: '', order: ''},
        { axis:'[hide]', name: '', ncofficients: '', lower_end: '', upper_end: '', order: ''},
        { axis:'[hide]', name: '', ncofficients: '', lower_end: '', upper_end: '', order: ''},
      ]
    }
  },
}
</script>

<style scoped>
.run-table-view {
  background: #F6F6F6;
  border-radius: 10px; 
  padding: 12px !important;
  height: 300px;
  overflow-y: scroll;
}

input[type=text] {
  width: 100%;
  background: #FFFFFF;
  border: 0.5px solid #838383;
  box-sizing: border-box;
  border-radius: 5px;
}

</style>