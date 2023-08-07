<template>
  <v-container fluid>
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
    <div id="clusterPlots">
      <!-- clusterPlots -->
      <v-row align="center">
        <v-col cols="3" align="center">
          <v-select
            class="pb-0"
            :items="barchartSortOptions"
            v-model="barchartSort_dropdown"
            label="sorting mode"
            dense
            hint=" sorting mode - chart/plot"
            persistent-hint
            single-line
            flat
          ></v-select>
        </v-col>
        <v-col cols="2" align="center">
          <v-switch
            class="ma-0"
            hide-details
            v-model="showNoise"
            label="with noise"
          ></v-switch>
        </v-col>
        <v-col cols="5" align="center">
          <v-radio-group class="mt-0" v-model="coloring" row hide-details>
            <v-tooltip top>
              <template v-slot:activator="{ on, attrs }">
                <v-radio
                  label="ID"
                  value="id"
                  v-bind="attrs"
                  v-on="on"
                ></v-radio>
              </template>
              <span>per Ligand</span>
            </v-tooltip>
            <v-tooltip top>
              <template v-slot:activator="{ on, attrs }">
                <v-radio
                  label="Mwt"
                  value="mwt"
                  v-bind="attrs"
                  v-on="on"
                ></v-radio>
              </template>
              <span>Molecular Weight</span>
            </v-tooltip>
            <v-tooltip top>
              <template v-slot:activator="{ on, attrs }">
                <v-radio
                  label="logP"
                  value="lipo"
                  v-bind="attrs"
                  v-on="on"
                ></v-radio>
              </template>
              <span>Octanol-water partition coefficient</span>
            </v-tooltip>
            <v-tooltip top>
              <template v-slot:activator="{ on, attrs }">
                <v-radio
                  label="Fsp3"
                  value="sp3"
                  v-bind="attrs"
                  v-on="on"
                ></v-radio>
              </template>
              <span> Fraction of sp³</span>
            </v-tooltip>
          </v-radio-group>
        </v-col>
      </v-row>

      <v-row>
        <v-col cols="10" id="wrapper" class="pr-0">
          <div id="tooltipCluster"></div>
          <div id="assigned"></div>
        </v-col>

        <v-col cols="2" id="scroll-wrapper-ligNames" class="pr-0 pl-0">
          <div id="legend" class="col pl-0 pr-0"></div>
        </v-col>
      </v-row>
    </div>
    <!-- end clusterPlots -->
  </v-container>
</template>
<script>
import {
  drawWrapper,
  prepareIncomingData,
  generateStackingData,
} from "@/assets/js/cluster";
import { mapState } from "vuex";

export default {
  name: "clustersComponent",
  data: () => ({
    loading: false,
    coloring: "mwt",
    lgdEnergies: {},
    sortValue: "Sort by Score",
    //showNoise: false,
  }),
  computed: {
    ...mapState({
      clusterData_recv: (state) => state.clusterData_recv,
      sorting: (state) => state.sorting,
      isFullscreen: (state) => state.isFullscreen,
    }),
    barchartSortOptions: {
      get: function () {
        return [
          { text: "Mean Energy/Score", value: "score" },
          { text: "Number of Ligands", value: "ligands" },
          { text: "Cluster ID", value: "clusters" },
        ];
      },
    },
    barchartSort_dropdown: {
      get: function () {
        return "clusters";
      },
      set: function (val) {
        console.log("barchartSortOption:", val);
        this.$store.commit("setSorting", val);
      },
    },
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
  },

  // Setup a Listener for one or more variables of "data", i.e. slider1 references to data -> slider1
  watch: {
    clusterData_recv: function () {
      if (this.clusterData_recv != null) {
        prepareIncomingData(this.clusterData_recv);
        generateStackingData();
        drawWrapper("#assigned");
        // no need for calling  drawWrapper("#assigned");
        // this watcher leads to a change of the mapped variable "sorting"
        // this leads to a call of  drawWrapper("#assigned");
      }
    },
    showNoise: function () {
      //document.getElementById("sortButton").innerHTML = " Sort by Energy";
      //document.getElementById("sortButton").value = "s";
      this.$store.commit("setColorMapSVGs", {});
      drawWrapper("#assigned");
    },
    coloring: function () {
      console.log("WATCH COLORING", this.coloring);
      this.$store.commit("setStackedBar_coloring", this.coloring);
      drawWrapper("#assigned");
    },
    sorting: function () {
      //console.log("sorting watch");
      drawWrapper("#assigned");
    },
    isFullscreen: function () {
      console.log("isFullscreen watch");
      setTimeout(() => {
        drawWrapper("#assigned");
      }, 1000);
    },
  },

  mounted() {
    //console.log(this.clusterData_recv)
    if (!(this.clusterData_recv == null)) {
      drawWrapper("#assigned");
    }
  },

  // methods of this vue component access via this.<funcName>
  methods: {},
};
</script>

<style></style>
