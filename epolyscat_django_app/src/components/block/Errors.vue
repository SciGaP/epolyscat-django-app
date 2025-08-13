<template>
  <div class="w-100 p-2">
    <b-alert v-for="(error, errorIndex) in errors" :key="errorIndex" show dismissible :variant="error.variant">
      <div class="w-100">
        <strong>{{ error.title }}</strong>
        <small> {{ error.description }} </small>
        <small v-if="error.source && error.source.stack" v-b-toggle="`error-stack-${errorIndex}`">
          <b-button size="sm" variant="link" class="when-closed">
            <small>Show more.</small>
          </b-button>
        </small>
      </div>
      <div class="w-100" v-if="error.source && error.source.stack">
        <b-collapse :id="`error-stack-${errorIndex}`" accordion="error-stack">
          <div class="overflow-auto">
            <small style="white-space: pre;">{{ error.source.stack }}</small>
          </div>
          <small v-if="error.source && error.source.stack" v-b-toggle="`error-stack-${errorIndex}`">
            <b-button size="sm" variant="link" class="when-open">
              <small>Show less.</small>
            </b-button>
          </small>
        </b-collapse>
      </div>
    </b-alert>
  </div>
</template>

<script>
export default {
  name: "Errors",
  props: {
    errors: {
      // {title, description, variant, source}
      default() {
        return [];
      }
    }
  }
}
</script>

<style scoped>
.collapsed > .when-open,
.not-collapsed > .when-closed {
  display: none;
}
</style>