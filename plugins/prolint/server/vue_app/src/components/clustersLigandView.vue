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
    <!--
    <v-navigation-drawer
      class="rounded"
      v-model="internalDrawer"
      absolute
      bottom
      temporary
      width="25%"
    >
      <component :is="navDrawer"></component>
    </v-navigation-drawer>
    -->

    <v-row>
      <v-col
        cols="6"
        class="d-flex pa-2 pr-1 flex-column"
        id="dockingOverviewComponent_all"
      >
        <!-- was col = 8-->
        <v-card class="mx-auto my-auto colCard" id="dockingOverviewComponent">
          <div id="LargeDiv">
            <div id="aboveFunctionalGroups">
              <!--- Larger Div --->
              <div class="pa-2">
                <v-row class="pt-2 pb-2">
                  <v-btn-toggle
                    v-model="toggle_view"
                    color="primary darken-3"
                    rounded
                    group
                  >
                    <v-btn @click="current = 'clustersComponent'">
                      Cluster
                    </v-btn>
                    <!-- <v-btn @click="current = 'ligandComponent'"> Ligand </v-btn>-->
                    <v-btn @click="current = 'heatmapComponent'">
                      Heatmap
                    </v-btn>
                  </v-btn-toggle>
                  <div id="infoCircle" class="pt-2">
                    <infoOverlayComponent
                      v-if="current == 'clustersComponent'"
                      class="pt-1"
                      header="Docking Overview:"
                    >
                    </infoOverlayComponent>
                    <infoOverlayComponent
                      v-if="current == 'ligandComponent'"
                      class="pt-1"
                      header="Docking Overview - LigandBars:"
                    >
                    </infoOverlayComponent>
                    <infoOverlayComponent
                      v-if="current == 'heatmapComponent'"
                      class="pt-1"
                      header="Docking Overview - Heatmap:"
                    >
                    </infoOverlayComponent>
                  </div>
                </v-row>
              </div>
              <!--<v-btn type="button" id="switchButton" v-model="clusterView" @click="switchView">Switch Cluster/Ligand View</v-btn>-->
              <v-card class="pa-0">
                <component :is="current"></component>
              </v-card>
            </div>
            <cluster-func-groups-component
              id="Functional_Groups"
              v-if="current == 'clustersComponent'"
              class="pb-0 pl-0 pr-0"
            ></cluster-func-groups-component>
          </div>
          <!-- end LargeDiv -->
        </v-card>
      </v-col>

      <!--- Medium Div --->
      <v-col cols="6" class="d-flex pa-2 pl-1 flex-column">
        <v-row>
          <v-col cols="12" class="d-flex pb-1 flex-column" id="clusterViewID">
            <!-- was col = 8-->
            <v-card class="mx-auto my-auto pa-0 pb-0 colCard">
              <v-card-title class="pb-0">
                Cluster Statistics:
                {{
                  clusterID == "All"
                    ? "All Clusters with Noise Poses"
                    : clusterID == "All_withoutNoise"
                    ? "All Clusters"
                    : "Cluster " + clusterID
                }}
                <infoOverlayComponent header="Statistics View:">
                </infoOverlayComponent>
                <v-spacer></v-spacer>
                <p class="caption">
                  filter mode:
                  <strong
                    :class="[
                      checkB_fgsClust == true
                        ? 'red--text text--lighten-2'
                        : 'blue--text text--lighten-2',
                    ]"
                  >
                    {{
                      checkB_fgsClust == true
                        ? "Functional Group Cluster"
                        : "Manual Treeview Selection"
                    }}
                  </strong>
                </p>
              </v-card-title>

              <component :is="'tooltipPLIPComponent'"></component>
              <component :is="'tooltipSkeletalFormulaComponent'"></component>

              <!--
              <div>

                <v-btn-toggle
                    v-model="toggle_table"
                    tile
                    color="primary darken-3"
                    group
                >
                  <v-btn @click="currentTable = 'clusterStatsComponent'">
                    Cluster Statistics
                  </v-btn>



                  <v-btn @click="currentTable = 'clusterFuncGroupsComponent'">
                    Cluster Functional Groups
                  </v-btn>


                </v-btn-toggle>

                <component :is="currentTable"></component>


              </div>
              -->

              <clusterStatsComponent class="pt-0"></clusterStatsComponent>
            </v-card>
          </v-col>
        </v-row>
        <v-row>
          <v-col cols="12" class="d-flex pt-0 pb-3 flex-column">
            <v-card class="mx-auto my-auto colCard">
              <v-card-title class="pb-1">
                Ligand: {{ ligandTableTitle }}
                <infoOverlayComponent header="Ligand View:">
                </infoOverlayComponent>
                <v-spacer></v-spacer>
                <!--<p class="caption">
                  color mode:
                  <strong
                    :class="[
                      checkB_fgsClust == true
                        ? 'red--text text--lighten-2'
                        : 'blue--text text--lighten-2',
                    ]"
                  >
                    {{
                      checkB_fgsClust == true
                        ? " Functional Group Cluster"
                        : "Functional Groups View"
                    }}
                  </strong>
                </p>-->
              </v-card-title>

              <v-row>
                <v-col cols="7" class="d-flex flex-column">
                  <v-card class="mx-auto my-auto colCard">
                    <v-card-title
                      class="pt-0"
                      style="
                        font-size: 11pt;
                        color: #7d7a7a;
                        font-weight: normal;
                      "
                      >Binding Modes

                      <v-spacer></v-spacer>
                      <v-tooltip right content-class="tooltipTable">
                        <template v-slot:activator="{ on, attrs }">
                          <v-btn
                            class="ma-1 pa-1"
                            outlined
                            small
                            color="indigo"
                            v-bind="attrs"
                            v-on="on"
                            id="colorLegendForces"
                          >
                            <div class="caption">
                              Interaction Types

                              <v-icon outlined small color="blue darken-2">
                                mdi-message-text dsd
                              </v-icon>
                            </div>
                          </v-btn>
                        </template>
                        <div>
                          <table>
                            <tr>
                              <td>
                                <v-icon
                                  class="pa-0 mx-0"
                                  dark
                                  left
                                  style="color: rgb(103.785, 49.98, 167.79)"
                                >
                                  mdi-circle
                                </v-icon>
                              </td>
                              <td>
                                <a style="color: black"> hydrogen bond</a>
                              </td>
                            </tr>
                            <tr>
                              <td>
                                <v-icon
                                  class="pa-0 mx-0"
                                  dark
                                  left
                                  style="color: rgb(55.08, 125.97, 184.11)"
                                >
                                  mdi-circle
                                </v-icon>
                              </td>
                              <td>
                                <a style="color: black"> halogen bond</a>
                              </td>
                            </tr>
                            <tr>
                              <td>
                                <v-icon
                                  class="pa-0 mx-0"
                                  dark
                                  left
                                  style="color: rgb(77.01, 174.93, 73.95)"
                                >
                                  mdi-circle
                                </v-icon>
                              </td>
                              <td>
                                <a style="color: black">
                                  hydrophobic interaction</a
                                >
                              </td>
                            </tr>
                            <tr>
                              <td>
                                <v-icon
                                  class="pa-0 mx-0"
                                  dark
                                  left
                                  style="color: rgb(166, 86, 40)"
                                >
                                  mdi-circle
                                </v-icon>
                              </td>
                              <td>
                                <a style="color: black"> metal complex</a>
                              </td>
                            </tr>
                            <tr>
                              <td>
                                <v-icon
                                  class="pa-0 mx-0"
                                  dark
                                  left
                                  style="color: rgb(255, 126.99, 0)"
                                >
                                  mdi-circle
                                </v-icon>
                              </td>
                              <td>
                                <a style="color: black">
                                  &pi;-cation interaction</a
                                >
                              </td>
                            </tr>
                            <tr>
                              <td>
                                <v-icon
                                  class="pa-0 mx-0"
                                  dark
                                  left
                                  style="color: rgb(245, 221, 25)"
                                >
                                  mdi-circle
                                </v-icon>
                              </td>
                              <td>
                                <a style="color: black">&pi;-stack</a>
                              </td>
                            </tr>
                            <tr>
                              <td>
                                <v-icon
                                  class="pa-0 mx-0"
                                  dark
                                  left
                                  style="color: rgb(226.95, 26.01, 28.05)"
                                >
                                  mdi-circle
                                </v-icon>
                              </td>
                              <td>
                                <a style="color: black"> salt bridge</a>
                              </td>
                            </tr>
                          </table>
                        </div>
                      </v-tooltip>
                    </v-card-title>

                    <div>
                      <!--
                      <v-card-title>
                        <v-text-field
                            v-model="tableSearch2"
                            append-icon="mdi-magnify"
                            label="Search Table"
                            single-line
                            hide-details
                            dense
                        ></v-text-field>
                      </v-card-title>
-->

                      <v-data-table
                        dense
                        :headers="tableHeaders2"
                        :items="tableData2"
                        :items-per-page="-1"
                        :search="tableSearch2"
                        :item-class="itemRowBackground"
                        @click:row="updateStateIDs"
                        id="lgdTable"
                        hide-default-footer
                        fixed-header
                        no-data-text="no ligand selected"
                        class="scrollableStyle_modelTable"
                      >
                        <!-- adds tooltip to table headers -->
                        <template
                          v-for="h in tableHeaders2"
                          v-slot:[`header.${h.value}`]="{ header }"
                        >
                          <v-tooltip top :key="h.value">
                            <template v-slot:activator="{ on }">
                              <span v-on="on">{{ header.text }}</span>
                            </template>
                            <span>{{ header.info }}</span>
                          </v-tooltip>
                        </template>

                        <template v-slot:[`item.mdl_id`]="{ item }">
                          <td>
                            <v-icon
                              class="pa-0 mx-0"
                              dark
                              left
                              dense
                              style="color: #9e9e9e"
                              @click="updateStateIDs_showFig(item.mdl_id)"
                            >
                              mdi-camera-outline
                            </v-icon>
                            {{ item.mdl_id }}
                          </td>
                        </template>

                        <template v-slot:[`item.hbond`]="{ item }">
                          <v-row>
                            <!-- HBonds -->
                            <v-icon
                              class="pa-0 mx-0"
                              v-if="item.hbond[0]"
                              dark
                              left
                              style="color: rgb(103.785, 49.98, 167.79)"
                            >
                              mdi-circle
                            </v-icon>
                            <v-icon
                              class="pa-0 mx-0"
                              v-else
                              dark
                              left
                              style="color: #f5f5f5"
                            >
                              mdi-circle
                            </v-icon>

                            <!-- halogenBonds -->
                            <v-icon
                              class="pa-0 mx-0"
                              v-if="item.hbond[1]"
                              dark
                              left
                              style="color: rgb(55.08, 125.97, 184.11)"
                            >
                              mdi-circle
                            </v-icon>
                            <v-icon
                              class="pa-0 mx-0"
                              v-else
                              dark
                              left
                              style="color: #f5f5f5"
                            >
                              mdi-circle
                            </v-icon>

                            <!-- hydrophobicInteractions -->
                            <v-icon
                              class="pa-0 mx-0"
                              v-if="item.hbond[2]"
                              dark
                              left
                              style="color: rgb(77.01, 174.93, 73.95)"
                            >
                              mdi-circle
                            </v-icon>
                            <v-icon
                              class="pa-0 mx-0"
                              v-else
                              dark
                              left
                              style="color: #f5f5f5"
                            >
                              mdi-circle
                            </v-icon>

                            <!-- metalComplexes -->
                            <v-icon
                              class="pa-0 mx-0"
                              v-if="item.hbond[3]"
                              dark
                              left
                              style="color: rgb(166, 86, 40)"
                            >
                              mdi-circle
                            </v-icon>
                            <v-icon
                              class="pa-0 mx-0"
                              v-else
                              dark
                              left
                              style="color: #f5f5f5"
                            >
                              mdi-circle
                            </v-icon>

                            <!-- piCationInteractions -->
                            <v-icon
                              class="pa-0 mx-0"
                              v-if="item.hbond[4]"
                              dark
                              left
                              style="color: rgb(255, 126.99, 0)"
                            >
                              mdi-circle
                            </v-icon>
                            <v-icon
                              class="pa-0 mx-0"
                              v-else
                              dark
                              left
                              style="color: #f5f5f5"
                            >
                              mdi-circle
                            </v-icon>

                            <!-- piStacks -->
                            <v-icon
                              class="pa-0 mx-0"
                              v-if="item.hbond[5]"
                              dark
                              left
                              style="color: rgb(245, 221, 25)"
                            >
                              mdi-circle
                            </v-icon>
                            <v-icon
                              class="pa-0 mx-0"
                              v-else
                              dark
                              left
                              style="color: #f5f5f5"
                            >
                              mdi-circle
                            </v-icon>

                            <!-- saltBridges -->
                            <v-icon
                              class="pa-0 mx-0"
                              v-if="item.hbond[6]"
                              dark
                              left
                              style="color: rgb(226.95, 26.01, 28.05)"
                            >
                              mdi-circle
                            </v-icon>
                            <v-icon
                              class="pa-0 mx-0"
                              v-else
                              dark
                              left
                              style="color: #f5f5f5"
                            >
                              mdi-circle
                            </v-icon>
                          </v-row>
                        </template>
                        <!--
                        <template v-slot:[`item.hbond`]="{ item }">
                          <v-icon v-if="item.hbond" dark left color="green">
                            mdi-check
                          </v-icon>
                          <v-icon v-else dark left color="red">
                            mdi-close
                          </v-icon>
                        </template>

                        <template v-slot:[`item.saltBridges`]="{ item }">
                          <v-icon v-if="item.saltBridges" dark left color="green">
                            mdi-check
                          </v-icon>
                          <v-icon v-else dark left color="red">
                            mdi-close
                          </v-icon>
                        </template>

                       <template v-slot:[`item.hydrophobicInteractions`]="{ item }">
                          <v-icon v-if="item.hydrophobicInteractions" dark left color="green">
                            mdi-check
                          </v-icon>
                          <v-icon v-else dark left color="red">
                            mdi-close
                          </v-icon>
                        </template>

                        <template v-slot:[`item.piCationInteractions`]="{ item }">
                          <v-icon v-if="item.piCationInteractions" dark left color="green">
                            mdi-check
                          </v-icon>
                          <v-icon v-else dark left color="red">
                            mdi-close
                          </v-icon>
                        </template>

                        <template v-slot:[`item.piStacks`]="{ item }">
                          <v-icon v-if="item.piStacks" dark left color="green">
                            mdi-check
                          </v-icon>
                          <v-icon v-else dark left color="red">
                            mdi-close
                          </v-icon>
                        </template>

                        <template v-slot:[`item.halogenBonds`]="{ item }">
                          <v-icon v-if="item.halogenBonds" dark left color="green">
                            mdi-check
                          </v-icon>
                          <v-icon v-else dark left color="red">
                            mdi-close
                          </v-icon>
                        </template>

                        <template v-slot:[`item.metalComplexes`]="{ item }">
                          <v-icon v-if="item.metalComplexes" dark left color="green">
                            mdi-check
                          </v-icon>
                          <v-icon v-else dark left color="red">
                            mdi-close
                          </v-icon>
                        </template>
-->
                      </v-data-table>
                    </div>
                  </v-card>
                </v-col>

                <v-col cols="5" class="d-flex flex-column">
                  <v-tooltip top>
                    <template v-slot:activator="{ on, attrs }">
                      <v-card class="mx-auto my-auto colCard">
                        <v-card-subtitle>Functional Groups </v-card-subtitle>
                        <div id="scrollableStyle_ligandTreeView">
                          <v-treeview
                            dense
                            hoverable
                            :active.sync="active"
                            activatable
                            :items="treeviewItems"
                            v-bind="attrs"
                            v-on="on"
                          >
                            <template v-slot:append="{ item }">
                              <div>
                                {{ item.amount }}
                              </div>
                            </template>
                          </v-treeview>
                        </div>
                      </v-card>
                    </template>
                  </v-tooltip>
                </v-col>
              </v-row>
            </v-card>
          </v-col>
        </v-row>
      </v-col>
    </v-row>
  </v-container>
</template>

<script>
import infoOverlayComponent from "./infoOverlayComponent";
import tooltipPLIPComponent from "./tooltipPLIPComponent";
import tooltipSkeletalFormulaComponent from "./tooltipSkeletalFormulaComponent";
import ligandComponent from "./ligandComponent";
import clustersComponent from "./clustersComponent";
import clusterStatsComponent from "./clusterStatsComponent";
import clusterFuncGroupsComponent from "./clusterFuncGroupsComponent";
import heatmapComponent from "./heatmapComponent";
//import navDrawer from './navDrawer';
import { mapState, mapGetters } from "vuex";
import fgClassifications from "@/assets/js/fgClassifications";
import {
  updateTable,
  highlightCluster,
  updateLgdTableData,
} from "../assets/js/cluster";
import { store } from "../main";
import { scrollParentToChild } from "./componentUtils";

//import {functionalGroupChart} from '@/assets/js/functionalGroupChart'
export default {
  // name of the component
  name: "clustersLigandView",
  props: ["drawer"],
  components: {
    ligandComponent,
    clustersComponent,
    heatmapComponent,
    clusterStatsComponent,
    clusterFuncGroupsComponent,
    tooltipPLIPComponent,
    infoOverlayComponent,
    tooltipSkeletalFormulaComponent,
    //navDrawer,
  },
  // data section of the Vue component. Access via this.<varName> .
  data: () => ({
    active: [],
    treeviewSearch: null,
    tooltipVar: false,
    current: "clustersComponent",
    currentTable: "clusterStatsComponent",
    loading: false,
    toggle_view: 0,
    toggle_table: 0,
    sorting: null,
    //sampleTextInput: "SAMPLE",
    //slider1: {label: "Slider1", val: [0,100], min: 0, max: 100},
    //radioTest: "radio1",
    //spinner1: 60,
    //spinner2: 2,
    tableSearch: "",
    tableSearch2: "",
    //ligandTableTitle: "Ligand",
    show_noise: false,
    functionalGroupChartObj: null,
    navDrawer: "navDrawer",
    infoBoxSegement: false,
  }),
  computed: {
    ...mapState({
      tableData: (state) => state.tableData,
      tableHeaders: (state) => state.tableHeaders,
      tableData2: (state) => state.tableData2,
      tableHeaders2: (state) => state.tableHeaders2,
      submitVal: (state) => state.submitVal,
      lgdEnergies: (state) => state.lgdEnergies,
      maxEnergy: (state) => state.maxEnergy,
      minEnergy: (state) => state.minEnergy,
      gblMdlID: (state) => state.gblMdlID,
      scrollToRow: (state) => state.modelID,
      ligandID: (state) => state.ligandID,
      modelID: (state) => state.modelID,
      clusterID: (state) => state.clusterID,
      colorScale: (state) => state.colorScale,
      ligandFunctionalGroups: (state) => state.funcGroupData_recv,
      selectedFuncGroup: (state) => state.selectedFuncGroup,
      largeLigand: (state) => state.largeLigand,
      currLigandFuncDat: (state) => state.currLigandFuncDat,
      functionalGroupsIDsToWord: (state) => state.functionalGroupsIDsToWord,
      mdlsInSelClust: (state) => state.mdlsInSelectedClust,
      ligandTableTitle: (state) => state.lgdTableTitle,
      curSelected_FGS_clusterData: (state) => state.curSelected_FGS_clusterData,
      checkB_fgsClust: (state) => state.checkB_fgsClust,
      supressSendData: (state) => state.supressSendData,
      clusterData_prepared: (state) => state.clusterData_prepared,
      selectionDataHash: (state) => state.selectionDataHash,
    }),
    ...mapGetters(["colorScale", "getZincSVG", "getFunctionalGroupPng"]),
    treeviewItems: {
      //TODO: Specificaiton tells us that objects are unordered atm this func relies on order
      get: function () {
        let outArr = [];
        if (this.currLigandFuncDat != null) {
          //console.log("TreeviewPart")

          let id = 1;
          let funcDatCopy = JSON.parse(JSON.stringify(this.currLigandFuncDat));
          //console.log("FuncGroups: ", funcDatCopy)
          let funcDatCopyKeys = Object.keys(funcDatCopy).reverse();
          while (Object.values(funcDatCopy).some((e) => e > 0)) {
            funcDatCopyKeys.forEach((element) => {
              while (funcDatCopy[element] > 0) {
                let tarArray = [];
                this.functionalGroupRecursion(tarArray, element);
                //console.log(tarArray)
                tarArray.forEach((elementInner) => {
                  funcDatCopy[elementInner] = funcDatCopy[elementInner] - 1;
                });
                //console.log(funcDatCopy)
                outArr.push(this.formatTreeviewEntry(tarArray, id));
                id += tarArray.length;
              }
            });
          }
          outArr = this.mergeDuplicates(outArr);
          outArr = this.recursiveMerge(outArr);
          //console.log(outArr)
        }
        return outArr;
      },
    },
    internalDrawer: {
      get: function () {
        return this.drawer;
      },
      set: function (newInternalDrawer) {
        this.$emit("update:drawer", newInternalDrawer);
      },
    },
  },
  // Setup a Listener for one or more variables of "data", i.e. slider1 references to data -> slider1
  watch: {
    scrollToRow: function () {
      //https://www.appsloveworld.com/vuejs/100/76/scroll-to-a-programmatically-selected-row-in-a-v-data-table
      /** @noinspection */
      setTimeout(() => {
        const row = document.getElementsByClassName(
          "highlightRow_selectedModelID"
        )[0];
        if (row) {
          const parent = document.getElementById("lgdTable").firstChild;
          scrollParentToChild(parent, row);
        }
      }, 10);
    },
    //spinner1: function(){console.log(this.spinner1)},
    //spinner2: function(){console.log(this.spinner2)},
    gblMdlID: function () {
      this.sendGblMdlID(this.gblMdlID);
    },
    ligandID: function () {
      store.commit("setLgdTableTitle", this.ligandID);
      //store.commit("setLargeLigand", store.getters.getZincSVG(preparedData.lgdNameID[ligandID]))
      store.commit("setCurrLigandFuncDat", this.ligandID);
      updateLgdTableData();
    },
    selectionDataHash: function () {
      if (!this.supressSendData) {
        if (isNaN(this.clusterID)) {
          this.sendClusterLigandModelID(-1, this.ligandID, this.modelID);
        } else {
          this.sendClusterLigandModelID(
            this.clusterID,
            this.ligandID,
            this.modelID
          );
        }
      }
    },
    ligandModelID: function () {},
    clusterID: function () {
      highlightCluster();
      updateTable();
    },
    clusterData_prepared: function () {
      if (this.clusterData_prepared != null) {
        updateTable();
      }
    },
    active: function () {
      this.searchTreeviewEntry(this.treeviewItems, this.active[0]);
      let val = this.treeviewSearch;
      this.setCurrentFuncGroup(val);
    },
    curSelected_FGS_clusterData: function () {
      updateLgdTableData();
    },
    checkB_fgsClust: function () {
      updateLgdTableData();
    },
  },

  // functions to call on mount (after DOM etc. is built)
  mounted() {
    //this.$vuetify.theme.themes.light.primary = '#2196F3';
  },

  // methods of this vue component access via this.<funcName>
  methods: {
    scrollParentToChild(parent, child) {
      // Where is the parent on page
      var parentRect = parent.getBoundingClientRect();
      // What can you see?
      var parentViewableArea = {
        height: parent.clientHeight,
        width: parent.clientWidth,
      };

      // Where is the child
      var childRect = child.getBoundingClientRect();
      // marco: "dif" as adjustment of scrolling sensitivity and length
      let dif = childRect.top - childRect.bottom;
      // Is the child viewable?
      var isViewable =
        childRect.top >= parentRect.top - dif * 2 &&
        childRect.bottom <= parentRect.top + parentViewableArea.height;
      // if you can't see the child try to scroll parent
      if (!isViewable) {
        // scroll by offset relative to parent
        parent.scrollTop =
          childRect.top +
          parent.scrollTop -
          parentRect.top -
          childRect.height +
          dif;
      }
    },
    searchTreeviewEntry(tar, query) {
      //console.log("TESTTAR", tar)
      //console.log("TESTQUERY", query)
      tar.forEach((elem) => {
        //console.log(elem)
        if (elem.id == query) {
          this.treeviewSearch = elem.val;
          return;
        } else if (elem.children) {
          this.searchTreeviewEntry(elem.children, query);
          return;
        }
      });
    },
    formatTreeviewEntry(unformattedArray, id) {
      let outObj = 0;
      for (let index = 0; index < unformattedArray.length; index++) {
        const element = unformattedArray[index];
        if (outObj == 0) {
          outObj = {
            id: id,
            val: element,
            name: this.functionalGroupsIDsToWord[element],
            amount: this.currLigandFuncDat[element],
            children: [],
          };
        } else {
          let tempObj = {
            id: id + unformattedArray.length - index,
            val: element,
            amount: this.currLigandFuncDat[element],
            name: this.functionalGroupsIDsToWord[element],
            children: [outObj],
          };
          outObj = tempObj;
        }
      }
      return outObj;
    },
    mergeDuplicates(originArray) {
      let outArray = [];
      let idBlacklist = [];
      for (let index = 0; index < originArray.length; index++) {
        let element = originArray[index];
        //console.log("MERGE elem", element)
        if (!idBlacklist.includes(element.id)) {
          for (
            let innerIndex = index + 1;
            innerIndex < originArray.length;
            innerIndex++
          ) {
            let compareElem = originArray[innerIndex];
            if (element.val == compareElem.val) {
              if (element.children) {
                element.children = element.children.concat(
                  compareElem.children
                );
              }
              idBlacklist.push(compareElem.id);
            }
          }

          outArray.push(element);
          //console.log("OUT ARR", element)
        }
      }
      //console.log("MERGE STEP",outArray)
      return outArray;
    },
    recursiveMerge(array) {
      array.forEach((element) => {
        if (element.children) {
          element.children = this.mergeDuplicates(element.children);
          //console.log("RECURSIVE MERGE STEP: ", element)
          if (element.children.length > 0) {
            this.recursiveMerge(element.children);
          }
        }
      }, array);
      return array;
    },
    functionalGroupRecursion(tar, query) {
      let val = fgClassifications[query];
      if (val == "root") {
        tar.push(query);
      } else {
        tar.push(query);
        this.functionalGroupRecursion(tar, val);
      }
    },
    itemRowBackground: function (item) {
      let rtv = "";
      if (item.mdl_id.toString() == (this.modelID + 1).toString()) {
        rtv += "highlightRow_selectedModelID ";
      }
      if (this.mdlsInSelClust[this.ligandID] != null) {
        for (let i = 0; i < this.mdlsInSelClust[this.ligandID].length; i++) {
          if (
            parseInt(item.mdl_id - 1) ===
            parseInt(this.mdlsInSelClust[this.ligandID][i])
          ) {
            if (parseInt(item.presModel_id) != -1) {
              return (
                rtv + "highlightRow_modelPresentInCluster_and_selectedFuncGroup"
              );
            } else {
              return rtv + "highlightRow_modelPresentInCluster";
            }
          }
        }
      }
      return rtv;
    },
    sleep(milliseconds) {
      return new Promise((resolve) => {
        setTimeout(() => {
          resolve();
        }, milliseconds);
      });
    },
    setCurrentFuncGroup(val) {
      this.$store.commit("setSelectedFuncGroup", parseInt(val));
    },
    sendGblMdlID(id) {
      this.$store.commit("setBlockUserInput", true);
      fetch("/sendMegaMolID?gblMdlID=" + id);

      //console.log("id-test",id);
    },
    sendClusterLigandModelID(clusterID, ligandID, modelID) {
      this.$store.commit("setBlockUserInput", true);
      //console.log("id-test",id0," ",id1);
      fetch(
        "/sendMegaMolID?clusterLigandModelID= &" +
          "clusterID=" +
          clusterID +
          "&" +
          "ligandID=" +
          ligandID +
          "&" +
          "modelID=" +
          modelID
      );
    },
    updateStateIDs(row) {
      let mdlID = parseInt(row["mdl_id"]);
      this.$store.commit("setModelID", mdlID - 1);

      //this.$store.commit('setClusterID', clustID);
      //console.log("Lgd Table Click set lgdID, mdlID, clusterID: ", this.ligandID, this.modelID, this.clusterID);
    },
    updateStateIDs_showFig(mdl_id) {
      let mdlID = parseInt(mdl_id);
      this.$store.commit("setModelID", mdlID - 1);
      let val = [];
      val.push(this.ligandTableTitle.split(" ")[1] + "_" + mdlID.toString());
      this.$store.commit("getProLigInterPNG", val);
      //this.$store.commit('setClusterID', clustID);
      //console.log("Lgd Table Click set lgdID, mdlID, clusterID: ", this.ligandID, this.modelID, this.clusterID);
    },
  },
};
</script>
