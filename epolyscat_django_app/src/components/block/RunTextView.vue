<template>
  <div class="run-text-view">
    <textarea class="run-text-class" name="run-text" id="" cols="30" rows="10" v-model="textValue">
    </textarea>
  </div>
</template>
<script>
export default {
  name: 'RunTextView',
  data() {
    return {
      textValue: 'The following input files describe specific applications: \n \nTRACER.P8C - Offline Infiltration Basin, Detention Pond in Series\nESM_L.P8C - Emerald Square Mall, Lower Watershed 2 Detention Ponds, Swale, 3 Wetland Cells in Series\nESM_U.P8C - Emerald Square Mall, Upper Watershed, Detention Pond, 3 Wetland Cells in Series\nHUNT.P8C - Hunt-Potowomut Watershed, RI. Watershed-Scale Application Illustrates use of Aquifer \nSample Precipitation Files\n\nPROVIDENCE_6987.PCP - Providence Airport, 1969-1987\nCITY_Y1Y2.PCP  - Precipitation Files for Other Cities between Y1 and Y2\n',
    }
  },
}
</script>
<style scoped>
.run-text-view {
  background: #F6F6F6;
  border-radius: 10px; 
  padding: 12px !important;
  height: 300px;
  overflow-y: scroll;
}

.run-text-class {
  width: 100%;
  background: #F6F6F6;
  border: none;
  resize: none;
}
</style>