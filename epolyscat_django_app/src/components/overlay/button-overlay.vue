<template>
  <b-overlay :show="show" rounded spinner-small :spinner-variant="variant" class="d-inline-block">
    <slot/>
  </b-overlay>
</template>

<script>
export default {
  name: "button-overlay",
  props: {
    show: {default: false},
    variant: {default: "primary"}
  },
}
</script>

<style scoped>

</style>