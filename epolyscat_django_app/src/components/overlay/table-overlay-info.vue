<template>
  <b-overlay :show="!data" rounded="sm" opacity="0" spinner-type="none">
    <b-skeleton-table v-if="data === null" :rows="rows" :columns="columns"
                      :table-props="{ bordered: true }"></b-skeleton-table>
    <slot name="empty" v-else-if="data.length === 0">
      <div class="w-100 text-center" style="color: #203A43;padding-top: 100px;">
        <p style="max-width: 600px;display: inline-block;">
          <img src="../../assets/not-found-icon.svg" style="width: 70px;margin-bottom: 30px;"/><br/>
          {{ emptyLabel }}.
        </p>
      </div>
    </slot>
    <slot v-else/>
  </b-overlay>
</template>

<script>
export default {
  name: "table-overlay-info",
  props: {
    data: {
      default() {
        return null
      }
    },
    rows: {},
    columns: {},
    emptyLabel: {
      default: "No data to be shown"
    }
  },
}
</script>

<style scoped>

</style>