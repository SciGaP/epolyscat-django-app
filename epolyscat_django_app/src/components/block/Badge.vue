<template>
    <span :class="`badge badge-${badgeVariant(status)} h-auto`">
        {{ status.toUpperCase() }}
    </span>
</template>

<script>
export default {
    props: ["status"],
    methods: {
        badgeVariant(status) {
            status = status.toUpperCase();

            if (["UNSUBMITTED"].indexOf(status) != -1) {
                return "secondary";
            } else if (["COMPLETED", "COMPLETE"].indexOf(status) != -1) {
                return "success";
            } else if (["CANCELING", "CANCELED", "SUSPENDED"].indexOf(status) != -1) {
                return "warning";
            } else if (["NON_CRITICAL_FAIL", "FAILED"].indexOf(status) != -1) {
                return "danger";
            } else {
                return "info";
            }
        }
    }
}
</script>