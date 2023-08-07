<template>
  <v-btn icon @click="clickAction">
    <v-icon
      class="mx-auto my-auto pa-2"
      v-if="showInfoIcons"
      dense
      :color="iconColor"
    >
      mdi-information
    </v-icon>
  </v-btn>
</template>

<script>
import { mapState } from "vuex";

export default {
  name: "infoOverlayComponent",
  props: {
    header: String,
    text: String,
    iconColor: String,
  },

  data: () => ({
    loading: false,
  }),
  computed: {
    ...mapState({
      showInfoIcons: (state) => state.showInfoIcons,

      // DOKU
      LigandView: (state) => state.LigandView,
      ClusterView: (state) => state.ClusterView,
      DockingOverview: (state) => state.DockingOverview,
      DockingOverview_LigandBars: (state) => state.DockingOverview_LigandBars,
      DockingOverview_Heatmap: (state) => state.DockingOverview_Heatmap,
      Functional_Groups_TreeView: (state) => state.Functional_Groups_TreeView,
      threeD_VisControls: (state) => state.threeD_VisControls,
    }),
    awesome: {
      get: function () {
        if (this.curProLigInterPNG) {
          return true;
        } else {
          return false;
        }
      },
    },
  },
  methods: {
    // DOKU
    /*
    threeD_VisControls
    DockingOverview_LigandBars
    DockingOverview_Heatmap
    Functional_Groups_TreeView
    LigandView
    ClusterView
    DockingOverview
    */
    clickAction() {
      this.$store.commit("setInfoOverlayState", true);
      if (this.header === "Docking Overview:") {
        this.$store.commit("setInfoOverlayContent", [
          this.header,
          "DockingOverview",
        ]);
      } else if (this.header === "Docking Overview - LigandBars:") {
        this.$store.commit("setInfoOverlayContent", [
          this.header,
          "DockingOverview_LigandBars",
        ]);
      } else if (this.header === "Docking Overview - Heatmap:") {
        this.$store.commit("setInfoOverlayContent", [
          this.header,
          "DockingOverview_Heatmap",
        ]);
      } else if (this.header === "Functional Groups View:") {
        this.$store.commit("setInfoOverlayContent", [
          this.header,
          "Functional_Groups_TreeView",
        ]);
      } else if (this.header === "Sidebar Menu:") {
        this.$store.commit("setInfoOverlayContent", [
          this.header,
          "threeD_VisControls",
        ]);
      } else if (this.header === "Statistics View:") {
        this.$store.commit("setInfoOverlayContent", [
          this.header,
          "ClusterView",
        ]);
      } else if (this.header === "Ligand View:") {
        this.$store.commit("setInfoOverlayContent", [
          this.header,
          "LigandView",
        ]);
      } else {
        this.$store.commit("setInfoOverlayContent", [
          "Header not defined",
          "dummy",
        ]);
      }
    },
  },
};
</script>
