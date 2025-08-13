<template>
  <div>
    <div>
      <b-button v-b-toggle.sidebar-right variant="link">
        <span class="text-vertical">
          Post Processing 1
          <b-icon icon="chevron-left" scale="0.5"></b-icon>
        </span>
      </b-button>
      <b-sidebar id="sidebar-right" title="Post Processing" right shadow>
        <PostProcessing/>
      </b-sidebar>
    </div>
    <div>
      <b-button v-b-toggle.sidebar-right variant="link">
        <span class="text-vertical">
          Post Processing 2
          <b-icon icon="chevron-left" scale="0.5"></b-icon>
        </span>
      </b-button>
    </div>
  </div>
</template>

<script>
import PostProcessing from "@/components/block/PostProcessing";

export default {
  name: "AppRightNav",
  components: {PostProcessing}
}
</script>

<style scoped>
.text-vertical {
  writing-mode: vertical-rl;
  text-orientation: mixed;
}
</style>