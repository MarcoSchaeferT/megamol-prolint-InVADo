<template>
  <v-container fluid class="pa-0">
    <!-- overlay https://stackoverflow.com/questions/53246403/vuetify-add-loading-layer-overlay-on-click-event  -->
    <v-overlay :value="loading" z-index="10">
      <v-container fluid fill-height>
        <v-layout justify-center align-center>
          <v-progress-circular
            :size="75"
            :width="10"
            indeterminate
            color="primary"
          >
          </v-progress-circular>
        </v-layout>
      </v-container>
    </v-overlay>

    <div id="vueSwitch_showNoise">
      <v-row no-gutters class="pt-4">
        <v-col cols="2">
          <v-switch
            class="ma-0 mb-2"
            hide-details
            v-model="showNoise"
            label="with noise"
          ></v-switch>
        </v-col>
        <v-col cols="3">
          <v-switch
            class="ma-0 pb-0 mb-2"
            hide-details
            v-model="dynamicRangeMin"
            label="dynamic max. color"
          ></v-switch>
        </v-col>
        <v-col cols="6">
          <div id="colormap"></div>
        </v-col>
      </v-row>
    </div>
    <div id="heatmapTotal">
      <!-- heatmap -->
      <div id="heatmapAxis"></div>
      <div id="heatmap">
        <div id="tooltipHeatmap"></div>
      </div>
    </div>
    <!-- end heatmap -->
  </v-container>
</template>
<!-- Load color palettes -->
<script src="https://d3js.org/d3-scale-chromatic.v1.min.js"></script>
<script>
import { drawHeatmap } from "@/assets/js/heatmap";
import { mapState } from "vuex";

export default {
  name: "heatmapComponent",
  data: () => ({
    loading: false,
    lgdEnergies: {},
  }),
  computed: {
    ...mapState({
      tableData: (state) => state.tableData,
      tableHeaders: (state) => state.tableHeaders,
      tableData2: (state) => state.tableData2,
      tableHeaders2: (state) => state.tableHeaders2,
      submitVal: (state) => state.submitVal,
      clusterData_recv: (state) => state.clusterData_recv,
      sorting: (state) => state.sorting,
      colorScale: (state) => state.colorScale,
    }),
    showNoise: {
      get() {
        return this.$store.state.showNoise;
      },
      set(val) {
        this.$store.commit("setShowNoise", val);
        this.$store.commit("calcMinEnergy"),
          this.$store.commit("calcMaxEnergy");
      },
    },
    dynamicRangeMin: {
      get() {
        return this.$store.state.dynamicEnergyMin;
      },
      set(val) {
        this.$store.commit("setDynamicEnergyMin", val),
          this.$store.commit("calcMinEnergy"),
          this.$store.commit("calcMaxEnergy");
      },
    },
  },
  // Setup a Listener for one or more variables of "data", i.e. slider1 references to data -> slider1
  watch: {
    showNoise: function () {
      drawHeatmap(
        this.clusterData_recv,
        "#heatmap",
        "#heatmapAxis",
        this.showNoise
      );
    },
    dynamicRangeMin: function () {
      drawHeatmap(
        this.clusterData_recv,
        "#heatmap",
        "#heatmapAxis",
        this.showNoise
      );
    },
  },
  mounted() {
    drawHeatmap(
      this.clusterData_recv,
      "#heatmap",
      "#heatmapAxis",
      this.showNoise
    );
  },
  // methods of this vue component access via this.<funcName>
  methods: {},
};
</script>

<style>
@import "../assets/css/heatmapStyles.css";
</style>
