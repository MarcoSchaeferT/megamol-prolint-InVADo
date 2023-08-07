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

    <div>
      <v-card-title class="pt-1 pb-0">
        <v-select
          :items="surfaceColoring_dropdownModes"
          v-model="ligandTableSorting_dropdown"
          label="filter by interaction type"
          outlined
          ><template #selection="{ item }">
            <span :style="'color:' + item.color">&#11044;</span>&nbsp;{{
              item.text
            }}
          </template>
          <template #item="{ item }">
            <span :style="'color:' + item.color">&#11044;</span>&nbsp;{{
              item.text
            }}
          </template>
        </v-select>
        <v-text-field
          v-model="tableSearch"
          append-icon="mdi-magnify"
          label="Search Table"
          outlined
        ></v-text-field>

        <v-spacer></v-spacer>
        <v-tooltip left content-class="tooltipTable">
          <template v-slot:activator="{ on, attrs }">
            <v-btn
              class="ma-1 pa-1"
              outlined
              small
              color="indigo"
              v-bind="attrs"
              v-on="on"
              id="colorLegend"
            >
              <div class="caption">
                Legend

                <v-icon outlined small color="blue darken-2">
                  mdi-message-text dsd
                </v-icon>
              </div>
            </v-btn>
          </template>

          <div
            class="tableColorMapContainer"
            v-html="getColorLegendTooltip"
          ></div>
        </v-tooltip>
      </v-card-title>
      <v-data-table
        dense
        :headers="tableHeaders"
        :items="tableData"
        :items-per-page="-1"
        :search="tableSearchReal"
        :custom-filter="ligandFilter"
        hide-default-footer
        :loading="checkForData"
        fixed-header
        item-key="lgd_id"
        :value="highlight_row"
        class="elevation-1 pb-0 pt-0 scrollableStyle_clusterStatistics"
        id="clusterTable"
        @click:row="updateLgdTable"
      >
        <!-- adds tooltip to table headers -->
        <template
          v-for="h in tableHeaders"
          v-slot:[`header.${h.value}`]="{ header }"
        >
          <v-tooltip top :key="h.value">
            <template v-slot:activator="{ on }">
              <span v-on="on">{{ header.text }}</span>
            </template>
            <span>{{ header.info }}</span>
          </v-tooltip>
        </template>
        <template v-slot:item.lgd_name="{ item }">
          <a
            style="color: darkcyan"
            :href="
              'https://zinc15.docking.org/substances/' +
              item.lgd_name.replace(/-[0-9]*/i, '')
            "
            target="_blank"
          >
            {{ item.lgd_name }} <v-icon small dense>mdi-open-in-new</v-icon>
          </a>
        </template>

        <template v-slot:item.svg="{ item }">
          <v-tooltip left content-class="tooltipTable">
            <template v-slot:activator="{ on }">
              <!-- V-HTML not the most secure solution-->
              <div
                v-on="on"
                class="tableSVGcontainer"
                v-html="getZincSVG(item.svg)"
              ></div>
            </template>
            <div v-html="getZincSVG(item.svg)"></div>
          </v-tooltip>
        </template>

        <template v-slot:item.local_min="{ item }">
          <v-tooltip left content-class="tooltipTable">
            <template v-slot:activator="{ on }">
              <v-chip v-on="on" :color="colorScale(item.local_min)" dark>
                {{ item.local_min }}
              </v-chip>
            </template>
            <div
              class="tableColorMapContainer"
              v-html="colormapTooltip(item.lgd_id, item.local_min, 'min')"
            ></div>
          </v-tooltip>
        </template>
        <template v-slot:item.local_max="{ item }">
          <v-tooltip right content-class="tooltipTable">
            <template v-slot:activator="{ on }">
              <v-chip v-on="on" :color="colorScale(item.local_max)" dark>
                {{ item.local_max }}
              </v-chip>
            </template>
            <div
              class="tableColorMapContainer"
              v-html="colormapTooltip(item.lgd_id, item.local_max, 'max')"
            ></div>
          </v-tooltip>
        </template>
      </v-data-table>
    </div>
  </v-container>
</template>

<script>
import { mapState, mapGetters } from "vuex";
import { drawColormapTooltip } from "@/assets/js/heatmap";
import { store } from "../main";
import { scrollParentToChild } from "./componentUtils";

//
// import {functionalGroupChart} from '@/assets/js/functionalGroupChart'
export default {
  // name of the component
  name: "clusterStatsComponent",

  // data section of the Vue component. Access via this.<varName> .
  data: () => ({
    loading: false,
    sortingForce: "",
    show_noise: false,
    largeLigand: "",
    functionalGroupChartObj: null,
  }),
  computed: {
    ...mapState({
      tableData: (state) => state.tableData,
      tableHeaders: (state) => state.tableHeaders,
      lgdEnergies: (state) => state.lgdEnergies,
      maxEnergy: (state) => state.maxEnergy,
      minEnergy: (state) => state.minEnergy,
      gblMdlID: (state) => state.gblMdlID,
      ligandID: (state) => state.ligandID,
      modelID: (state) => state.modelID,
      clusterID: (state) => state.clusterID,
      colorScale: (state) => state.colorScale,
      ligandFunctionalGroups: (state) => state.funcGroupData_recv,
      selectedFuncGroup: (state) => state.selectedFuncGroup,
      currLigandFuncDat: (state) => state.currLigandFuncDat,
      presentClusters: (state) => state.presentClustData,
      tableSearch: (state) => state.tableSearch,
      clusterData_recv: (state) => state.clusterData_recv,
      colorMapSVGs: (state) => state.colorMapSVGs,
      scrollToRow: (state) => state.ligandID,
      surfaceColoring_dropdownModesStore: (state) =>
        state.surfaceColoring_dropdownModesStore,
    }),

    surfaceColoring_dropdownModes: {
      get: function () {
        // eslint-disable-next-line vue/no-side-effects-in-computed-properties
        let tmp = [...this.surfaceColoring_dropdownModesStore];
        tmp.push({
          text: "no filter",
          value: 1,
          color: "rgb(255, 255, 255)",
          varName: "",
        });

        return tmp;
      },
    },
    ligandTableSorting_dropdown: {
      get: function () {
        return 1;
      },
      set: function (val) {
        for (let ele of this.surfaceColoring_dropdownModes) {
          if (ele.value == val) this.sortingForce = ele.varName;
        }
        if (val == 1) this.sortingForce = "";
        this.ligandFilter;
        console.log("sortingForce", this.sortingForce);
        //this.$store.commit("setTableSearch", "*");
      },
    },

    ...mapGetters(["colorScale", "getZincSVG", "getFunctionalGroupPng"]),
    highlight_row: {
      get() {
        let out = [{ lgd_id: this.ligandID.toString() }];
        return out;
      },
      set() {
        return 1;
      },
    },
    checkForData: {
      get() {
        if (this.tableData.length > 0) {
          return false;
        } else {
          return true;
        }
      },
    },
    tableSearch: {
      get() {
        return this.$store.state.tableSearch;
      },
      set(val) {
        this.$store.commit("setTableSearch", val);
      },
    },
    tableSearchReal: {
      get() {
        let searcher = "";
        if (
          this.$store.state.tableSearch == "" ||
          this.$store.state.tableSearch == null
        ) {
          searcher = this.sortingForce;
        } else {
          // searcher += this.sortingForce + ";";
          searcher += this.$store.state.tableSearch;
        }
        console.log("searcher", searcher);
        return searcher;
      },
    },
    currentlySelectedFunctionalGroup: {
      get: function () {
        if (this.selectedFuncGroup != null) {
          //console.log("COMPUTED TRIGGERED",this.selectedFuncGroup)
          let base64img = this.getFunctionalGroupPng(this.selectedFuncGroup);
          //console.log("image recieved", base64img)
          //let convertedImg = new Image()
          //convertedImg.src= base64img
          return base64img;
        }
        return null;
      },
    },
    getColorLegendTooltip() {
      return this.colorMapSVGs["dummy"];
    },
  },

  // functions to call on mount (after DOM etc. is built)
  mounted() {
    //this.$vuetify.theme.themes.light.primary = '#2196F3';
  },
  watch: {
    scrollToRow: function () {
      //https://www.appsloveworld.com/vuejs/100/76/scroll-to-a-programmatically-selected-row-in-a-v-data-table
      /** @noinspection */
      setTimeout(() => {
        const row = document.getElementsByClassName(
          "v-data-table__selected"
        )[0];
        const parent = document.getElementById("clusterTable").firstChild;
        scrollParentToChild(parent, row);
      }, 10);
    },
  },
  // methods of this vue component access via this.<funcName>
  methods: {
    ligandFilter(value, search, item) {
      if (search.includes(",")) {
        let ligandArray = search.split(",");
        if (value === this.sortingForce || this.sortingForce == "") {
          //console.log(value, " === ", this.sortingForce, item.lgd_id);
          if (ligandArray.includes(item.lgd_name) || ligandArray.length <= 1) {
            return value != null && search != null && typeof value === "string";
          }
        }
      } else {
        console.log("enter normal search");
        if (
          item.forceIsPresentLig.HBonds === this.sortingForce ||
          item.forceIsPresentLig.halogenBonds === this.sortingForce ||
          item.forceIsPresentLig.hydrophobicInteractions ===
            this.sortingForce ||
          item.forceIsPresentLig.metalComplexes === this.sortingForce ||
          item.forceIsPresentLig.piCationInteractions === this.sortingForce ||
          item.forceIsPresentLig.piStacks === this.sortingForce ||
          item.forceIsPresentLig.saltBridges === this.sortingForce ||
          this.sortingForce == ""
        ) {
          return (
            value != null &&
            search != null &&
            typeof value === "string" &&
            value.toString().indexOf(search) !== -1
          );
        } else {
          return false;
        }
      }
    },

    updateLgdTable(row) {
      //console.log("row",row);
      let lgdID = parseInt(row["lgd_id"]);
      let lgdName = row["lgd_name"];

      // Determine modelID with highest energy score
      // nrgs are saved as strings, convert to number first

      let nrgArr =
        this.presentClusters[lgdID]["clusters"][this.clusterID]["nrgs"].map(
          Number
        );
      let maxNrgModelIndex = nrgArr.indexOf(Math.min(...nrgArr));
      let maxNrgModelID =
        this.presentClusters[lgdID]["clusters"][this.clusterID]["models"][
          maxNrgModelIndex
        ];

      this.$store.commit("setLigandID", lgdID);
      this.$store.commit("setModelID", maxNrgModelID);

      this.$store.commit("setLgdTableTitle", lgdID);
      this.$store.commit("setLargeLigand", this.getZincSVG(lgdName));
      this.$store.commit("setCurrLigandFuncDat", lgdID);
      /*
      if(this.functionalGroupChartObj == null){
        this.functionalGroupChartObj = new functionalGroupChart(["#functionalGroupsHorizonalAxis","#functionalGroupsMain"], this.currLigandFuncDat)
      }else{
        this.functionalGroupChartObj.updateData(this.currLigandFuncDat)
      }
      */
      //console.log(lgdName)
      //console.log(this.ligandFunctionalGroups)
      //updateLgdTableData();
    },
    colormapTooltip(lgd_id, val, name) {
      let margin = {
        top: 20,
        bottom: 20,
        left: 60,
        right: 30,
      };
      if (this.colorMapSVGs["dummy"] == null) {
        let svgEntry = {};
        let margin = {
          top: 20,
          bottom: 20,
          left: 60,
          right: 30,
        };
        svgEntry["dummy"] = drawColormapTooltip(margin, 600, 0);
        this.$store.commit("setColorMapSVGs", svgEntry);
      }

      let entryName = this.clusterID.toString() + lgd_id.toString() + name;
      if (this.colorMapSVGs[entryName] == undefined) {
        let tmp = this.colorMapSVGs;
        tmp[entryName] = drawColormapTooltip(margin, 600, val);
        store.commit("setColorMapSVGs", tmp);
      }

      return this.colorMapSVGs[entryName];
    },
    sortingSelection(val) {
      console.log("current sorting selection: ", val);
    },
  },
};
</script>
