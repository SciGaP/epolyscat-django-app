export default {
  data() {
    return { page: 1, pageSize: 10, search: "", searchError: "" };
  },
  computed: {
    collectionQuery() {
      return { page: this.page, pageSize: this.pageSize, search: this.search.trim() };
    },
  },
  watch: {
    search() {
      this.page = 1;
      this[this.$options.collectionSelectionKey] = {};
    },
    page() {
      this[this.$options.collectionSelectionKey] = {};
    },
    pageSize() {
      this.page = 1;
    },
    collectionQuery() {
      this.refreshData();
    },
  },
  methods: {
    async refreshData() {
      const query = this.collectionQuery;
      this.searchError = "";
      try {
        await this.$store.dispatch(this.$options.collectionFetchAction, query);
      } catch (error) {
        if (JSON.stringify(query) === JSON.stringify(this.collectionQuery)) {
          this.searchError = "Unable to load this list. Please try again.";
        }
      }
    },
  },
  mounted() {
    this.refreshData();
  },
};
