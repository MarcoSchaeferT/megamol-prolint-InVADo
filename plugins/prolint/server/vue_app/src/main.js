import Vue from "vue";
import Vuex from "vuex";
import * as d3 from "d3";
import { functionalGroupsIDsToWord } from "@/assets/js/funcGroupLookUp";
import "@mdi/font/css/materialdesignicons.css";

import App from "./App.vue";
import vuetify from "./plugins/vuetify";
import { buildLgdInfoDict, parseFGSMapData } from "@/assets/js/cluster";
import { marked } from "marked";
import { rgb } from "d3";

function delay(time) {
  return new Promise((resolve) => setTimeout(resolve, time));
}

// the keys used to access the 'dataset' variable (the cluster- and ligand-data)
export const clusterDataKeys = {
  ligandCount: 0,
  modelCounts: 1,
  energies: 2,
  ligandIDs: 3,
  zincNames: 4,
  clusterAssignment: 5,
  clusterSizes: 6,
  // [7] ->  dataPath (path where svg, smi, checkmol, etc. are located)
  dataPath: 7,
  hbonds: 8,
  efficiency: 9,
  mwt: 10,
  logp: 11,
  fractioncsp3: 12,
  halogenBonds: 13,
  hydrophobicInteractions: 14,
  metalComplexes: 15,
  piCationInteractions: 16,
  piStacks: 17,
  saltBridges: 18,
};
const k = clusterDataKeys;
export const funcGroupDataKeys = {
  funcGroupsPerLigands: 0,
  funcGroupsPerClusters: 1,
};

Vue.config.productionTip = false;

Vue.use(Vuex);
Vue.use(vuetify, {
  iconfont: "mdi",
});
export const store = new Vuex.Store({
  state: {
    testMode: false,
    testModeRefreshTime: 60000,
    numberOf_FGStypeIDs: 204,
    clusterData_prepared: null,
    clusterData_recv: null,
    funcGroupData_recv: null,
    funcGroupData_inverted: null,
    functionalGroupsIDsToWord: functionalGroupsIDsToWord,
    tableHeaders: [
      {
        text: "ID",
        align: "start",
        sortable: true,
        value: "lgd_id",
        info: "ID",
      },
      {
        text: "Ligand Name",
        align: "start",
        sortable: true,
        value: "lgd_name",
        info: "ZINC name and web-link",
        width: "178px",
      },
      {
        text: "Structural Formula",
        align: "start",
        sortable: false,
        value: "svg",
        info: "structural formula",
      },
      {
        text: "Local Pose Count",
        align: "start",
        sortable: true,
        value: "mdl_cnt",
        info: "number of poses which are present in 'selected' cluster",
      },
      {
        text: "Local Min. Score",
        align: "start",
        sortable: true,
        value: "local_min",
        info: "min. score/max. energy for this ligand in 'selected' cluster",
      },
      {
        text: "Local Max. Score",
        align: "start",
        sortable: true,
        value: "local_max",
        info: "max. score/min. energy for this ligand in 'selected' cluster",
      },
      {
        text: "Average Pose Score",
        align: "start",
        sortable: true,
        value: "avg_mdl",
        info: "sum of all pose scores in total / total poses count (including noise poses)",
      },
      {
        text: "Present Clusters Count",
        align: "start",
        sortable: true,
        value: "clust_cnt",
        info: "number of clusters in which this ligand is present",
      },
      {
        text: "Current Cluster Average",
        align: "start",
        sortable: true,
        value: "curr_clust_avg",
        info: "sum of scores of poses which are present in 'selected' cluster / count of them",
      },
      {
        text: "Total Cluster Average",
        align: "start",
        sortable: true,
        value: "total_clust_avg",
        info: "sum of scores of poses which are present in 'any' found cluster / count of them",
      },
      {
        text: "HBonds",
        align: " d-none",
        sortable: true,
        value: "forceIsPresentLig.HBonds",
      },
      {
        text: "halogenBonds",
        align: " d-none",
        sortable: true,
        value: "forceIsPresentLig.halogenBonds",
      },
      {
        text: "hydrophobicInteractions",
        align: " d-none",
        sortable: true,
        value: "forceIsPresentLig.hydrophobicInteractions",
      },
      {
        text: "metalComplexes",
        align: " d-none",
        sortable: true,
        value: "forceIsPresentLig.metalComplexes",
      },
      {
        text: "piCationInteractions",
        align: " d-none",
        sortable: true,
        value: "forceIsPresentLig.piCationInteractions",
      },
      {
        text: "piStacks",
        align: " d-none",
        sortable: true,
        value: "forceIsPresentLig.piStacks",
      },
      {
        text: "saltBridges",
        align: " d-none",
        sortable: true,
        value: "forceIsPresentLig.saltBridges",
      },
    ],
    tableData: [],
    tableSearch: "",

    tableHeaders2: [
      {
        text: "Pose",
        align: "center",
        sortable: true,
        value: "mdl_id",
        width: "13%",
        info: "model/pose ID of the ligand",
      },
      {
        text: "Score",
        align: "start",
        sortable: true,
        value: "mdl_nrg",
        width: "1%",
        info: "score/free binding energy of ligand pose",
      },
      {
        text: "Ligand Efficiency",
        align: "start",
        sortable: true,
        value: "lgd_eff",
        width: "1%",
        info: "score/energy devided by number of non-hydrogen atoms",
      },
      {
        text: "Cluster-ID      ",
        align: "top",
        sortable: true,
        value: "presClust_id",
        width: "1%",
        info: "cluster/pocket ID the ligand binding pose belongs to",
      },
      {
        text: "Interaction Types",
        align: "start",
        sortable: true,
        value: "hbond",
        width: "35%",
        info: "interaction types of the ligand binding pose (hover INTERACTION TYPES button above on the right side )",
      },
    ],
    // {text: "salt bridge",align: "start", sortable: true ,value:"saltBridges"},
    // {text: "hyd. phob. Inter",align: "start", sortable: true ,value:"hydrophobicInteractions"},
    // {text: "pi+ inter",align: "start", sortable: true ,value:"piCationInteractions"},
    // {text: "pi inter",align: "start", sortable: true ,value:"piStacks"},
    //{text: "metal comp.",align: "start", sortable: true ,value:"metalComplexes"},
    //{text: "halogen bond",align: "start", sortable: true ,value:"halogenBonds"},],
    tableData2: [],

    stackingDat: null,
    currentTooltip_PLIP: "",

    functionalGroupsPNGs: {},
    submitVal: "0,0",
    nrgDict: null,
    lgdEnergies: "",
    maxEnergy: 0,
    minEnergy: 0,
    gblMdlID: 0,
    ligandID: 0,
    modelID: 0,
    clusterID: "All_withoutNoise",
    sorting: "clusters",
    ldgMdl: null,
    dynamicEnergyMin: true,
    showNoise: false,
    zincSVGs: {},
    proLigInterPNGs: {},
    selectedFuncGroup: null,
    largeLigand: null,
    currLigandFuncDat: null,
    clusterFuncData: null,
    mdlsInSelectedClust: {},
    lgdTableTitle: "Ligand",
    presentClustData: {},
    hbonds: {},
    halogenBonds: {},
    hydrophobicInteractions: {},
    metalComplexes: {},
    piCationInteractions: {},
    piStacks: {},
    saltBridges: {},

    efficiency: {},

    barScaleX: null,
    barScaleY: null,
    highlightBar: null,
    stackedBar_coloring: "mwt",

    // functional groups
    checkedFGSgroups: null,
    fgsMapData_recv: null,
    curSelected_FGS_clusterData: null,
    new_FGSclusterData: false,

    // treeView
    treeActives: [],

    // checkBoxes
    checkB_fgsClust: false,

    // pollFunctions Params
    refreshTime: 1000,

    // data hash
    GUIdataHash: -1,
    MegaMol_selectionDataHash: -1,

    supressSendData: false,

    // GUI data
    GUIfloatParams: [
      "modelClustering_eps",
      "modelClustering_minBindEnergy",
      "searchRadius_fgc",
      "clipPlaneDist",
    ],
    GUIIntParams: [
      "modelClustering_minPts",
      "minPTS_fgc",
      "starMenuCircleCnt",
      "starMenuSize",
      "starMenuTextSize",
      "subSphereSizefactor",
      "proteinColoringMode0",
      "proteinColoringMode1",
      "clipPlaneOrientation",
      "surfaceColoring_dropdown",
      "starMenuSorting_dropdown",
    ],
    GUIdata: [],
    surfaceColoring_dropdownModesStore: [
      {
        text: "hydrogen bonds",
        varName: "HBonds",
        value: 2,
        color: "rgb(103.785, 49.98, 167.79)",
      },
      {
        text: "halogen bonds",
        varName: "halogenBonds",
        value: 3,
        color: "rgb(55.08, 125.97, 184.11)",
      },
      {
        text: "hydrophobic interactions",
        varName: "hydrophobicInteractions",
        value: 4,
        color: "rgb(77.01, 174.93, 73.95)",
      },
      {
        text: "metal complexes",
        varName: "metalComplexes",
        value: 5,
        color: "rgb(166, 86, 40)",
      },
      {
        text: "π-cation interactions",
        varName: "piCationInteractions",
        value: 6,
        color: "rgb(255, 126.99, 0)",
      },
      {
        text: "π-stacks",
        varName: "piStacks",
        value: 7,
        color: "rgb(245, 221, 25)",
      },
      {
        text: "salt bridges",
        varName: "saltBridges",
        value: 8,
        color: "rgb(226.95, 26.01, 28.05)",
      },
    ],

    starMenuSortingForces_dropdownModesStore: [
      {
        text: "hydrogen bonds",
        varName: "HBonds",
        value: 2,
        color: "rgb(103.785, 49.98, 167.79)",
      },
      {
        text: "halogen bonds",
        varName: "halogenBonds",
        value: 3,
        color: "rgb(55.08, 125.97, 184.11)",
      },
      {
        text: "hydrophobic interactions",
        varName: "hydrophobicInteractions",
        value: 4,
        color: "rgb(77.01, 174.93, 73.95)",
      },
      {
        text: "metal complexes",
        varName: "metalComplexes",
        value: 5,
        color: "rgb(166, 86, 40)",
      },
      {
        text: "π-cation interactions",
        varName: "piCationInteractions",
        value: 6,
        color: "rgb(255, 126.99, 0)",
      },
      {
        text: "π-stacks",
        varName: "piStacks",
        value: 7,
        color: "rgb(245, 221, 25)",
      },
      {
        text: "salt bridges",
        varName: "saltBridges",
        value: 8,
        color: "rgb(226.95, 26.01, 28.05)",
      },
      {
        text: "best scores",
        varName: "bestScores",
        value: 9,
        color: "rgb(255, 255, 255)",
      },
    ],

    // GUIcolorScaleNames not in use --> just for docu reasons
    GUIcolorScaleNames: {
      interactionForce: "interactionForce",
      pocketCluster: "pocketCluster",
      funcGroupCluster: "funcGroupCluster",
    },
    GUIcolorScaleRanges: {
      interactionForce: [
        d3.rgb(0.4466 * 255, 0.6466 * 255, 0.9866 * 255, 1.0 * 255),
        d3.rgb(1.0 * 255, 0.02666 * 255, 1.0 * 255, 1.0 * 255),
      ],
      pocketCluster: [
        d3.rgb(1.0 * 255, 1.0 * 255, 1.0 * 255, 1.0 * 255),
        d3.rgb(0.1066 * 255, 1.0 * 255, 0.1681 * 255, 1.0 * 255),
      ],
      funcGroupCluster: [
        d3.rgb(1.0 * 255, 1.0 * 255, 1.0 * 255, 1.0 * 255),
        d3.rgb(1 * 255, 0.0 * 255, 0.0 * 255, 1.0 * 255),
      ],
    },

    GUIcolorScaleSVGs: {},

    showNavDrawer: false,
    showBlockUserInput: false,

    colorMapSVGs: {},

    selectionDataHash: 0,

    curProLigInterPNG: 0,
    curProLigInterPNG_name: "",
    showInfoIcons: true,
    showInfoOverlay: false,
    infoOverlayContent: ["Header", "text"],

    isFullscreen: false,

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
  },
  getters: {
    colorScale: (state) => {
      let scale;
      if (state.dynamicEnergyMin) {
        scale = d3
          .scaleSequential()
          .domain([-state.minEnergy, -state.maxEnergy])
          .interpolator(d3.interpolateInferno);
      } else {
        scale = d3
          .scaleSequential()
          .domain([0, -state.maxEnergy])
          .interpolator(d3.interpolateInferno);
      }
      return scale;
    },
    getZincSVG: (state) => (val) => {
      //const parser = new DOMParser();
      //console.log("Z-ID", val)
      //console.log("TEST-S",state.zincSVGs[val])
      //let svg = parser.parseFromString(state.zincSVGs[val],"image/svg+xml")
      //console.log("TEST-X",svg)
      return state.zincSVGs[val];
    },
    getFunctionalGroupPng: (state) => (val) => {
      return state.functionalGroupsPNGs[val];
    },
  },
  actions: {
    handleClusterData({ commit, state }) {
      //console.log("BEFORE FETCH: parsing cluster data");
      /* Datapackage
        [0]	->	ligandCount
        [1]	->	modelCount
        [2]	->	energies
        [3]	->	ligandIDs
        [4]	->	zincNames
        [5]	->	clusterAssignment
        [6]	->	clusterSizes
        [7] ->  dataPath (path where svg, smi, checkmol, etc. are located)
        [8] ->  hbonds
        [9] ->  efficiency
        [10] ->  mwt (molecular weight)
        [11] ->  logp (logP is a quantitative measure of lipophilicity)
        [12] ->  fractioncsp3 (parameter for drug-likeness)
        [13] ->  halogenBonds
        [14] ->  hydrophobicInteractions
        [15] ->  metalComplexes
        [16] ->  piCationInteractions
        [17] ->  piStacks
        [18] ->  saltBridges
      */

      return fetch("/getClusterData", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        //body: JSON.stringify({exportData: this.slider1.val})
      })
        .then((response) => response.json())
        .then((clusterData) => {
          return new Promise((resolve) => {
            console.log("AFTER FETCH: parsing cluster data");
            commit("setClusterData", clusterData.clusterData_recv);
            commit("calcMaxEnergy");
            commit("calcMinEnergy");
            commit("setSorting", state.sorting);
            commit(
              "setLgdEnergies",
              buildLgdInfoDict(clusterData.clusterData_recv, k.energies)
            );
            commit("setLgdMdl", clusterData.lgdMdl);
            commit(
              "setHbonds",
              buildLgdInfoDict(clusterData.clusterData_recv, k.hbonds)
            );
            commit(
              "setHalogenBonds",
              buildLgdInfoDict(clusterData.clusterData_recv, k.halogenBonds)
            );
            commit(
              "setHydrophobicInteractions",
              buildLgdInfoDict(
                clusterData.clusterData_recv,
                k.hydrophobicInteractions
              )
            );
            commit(
              "setMetalComplexes",
              buildLgdInfoDict(clusterData.clusterData_recv, k.metalComplexes)
            );
            commit(
              "setPiCationInteractions",
              buildLgdInfoDict(
                clusterData.clusterData_recv,
                k.piCationInteractions
              )
            );
            commit(
              "setPiStacks",
              buildLgdInfoDict(clusterData.clusterData_recv, k.piStacks)
            );
            commit(
              "setSaltBridges",
              buildLgdInfoDict(clusterData.clusterData_recv, k.saltBridges)
            );
            commit(
              "setEfficiency",
              buildLgdInfoDict(clusterData.clusterData_recv, k.efficiency)
            );
            resolve();
          });
        })
        .then(() => {
          //console.log("END: parsing cluster data")
          return fetch("/getOperationMode", {
            method: "POST",
            headers: {
              "Content-Type": "application/json",
            },
            //body: JSON.stringify({exportData: this.slider1.val})
          })
            .then((response) => response.json())
            .then((operationMode) => {
              return new Promise((resolve) => {
                let val = false;
                if (operationMode.operatioMode_recv === "testMode") {
                  val = true;
                }
                //console.log("testMode: ", val);
                commit("setOperationMode", val);
                if (val) commit("setPollTime", this.state.testModeRefreshTime);
                resolve();
              });
            });
        });
    },
    create_GUIcolorScaleSVGs({ commit }) {
      for (const [key, value] of Object.entries(
        this.state.GUIcolorScaleRanges
      )) {
        buildColorScaleSVG(value, key);
      }

      function buildColorScaleSVG(colorRange, scaleName) {
        let interpolator = d3.scaleLinear().domain([0, 1]).range(colorRange);
        const colorMapState = {
          rangeStart: 0,
          n: 200,
          colorRange,
          colorScale: d3.scaleSequential(interpolator),
        };
        let svgWidth = 325;
        let svgHeight = 13;
        let margin = {
          left: 10,
          right: 15,
          top: 10,
          bottom: 10,
        };
        let colorScaleInteractionForceSVG = d3
          .create("svg")
          .attr("width", svgWidth + margin.left + margin.right)
          .attr("height", svgHeight + margin.top + margin.bottom)
          .attr(
            "transform",
            "translate(" + margin.left + ",-" + margin.top + ")"
          );

        // draw color scale
        colorScaleInteractionForceSVG
          .selectAll()
          .data(d3.range(0, 1, 1 / colorMapState.n))
          .enter()
          .append("rect")
          .attr("class", "bars")
          .attr("y", 0)
          .attr("height", svgHeight)
          .attr("width", svgWidth / colorMapState.n + 1)
          .attr("x", (d, i) => margin.left + (i * svgWidth) / colorMapState.n)
          .style("fill", (d) => colorMapState.colorScale(d));
        //.style("fill", colorMapState.colorScale);

        // draw stroke around the color scale
        colorScaleInteractionForceSVG
          .append("rect")
          .attr("class", "bars")
          .attr("y", 0)
          .attr("height", svgHeight)
          .attr("width", svgWidth)
          .attr("x", margin.left + svgWidth / colorMapState.n - 1)
          .attr("stroke", "#000000")
          .attr("stroke-linecap", "butt")
          .attr("stroke-width", "1")
          .style("fill", rgb(1, 1, 1, 0));

        //x-axis scale
        let x = d3
          .scalePoint()
          .range([0, (svgWidth / colorMapState.n) * colorMapState.n])
          .domain(["min", "max"]);

        //create x-axis
        colorScaleInteractionForceSVG
          .append("g")
          .attr("class", "axis")
          .attr("transform", "translate(" + margin.left + "," + svgHeight + ")")
          .call(d3.axisBottom(x))
          .style("color", "#000000");

        var serializer = new XMLSerializer();
        var source = serializer.serializeToString(
          colorScaleInteractionForceSVG.node()
        );
        commit("setGUIcolorScaleSVGs", [scaleName, source]);
      }
    },

    pollFGSMapData({ commit, state }) {
      function f1() {
        let interval = setInterval(
          () => {
            return fetch("/getFGSMapData", {
              method: "POST",
              headers: {
                "Content-Type": "application/json",
              },
            })
              .then((response) => response.json())
              .then((FGSMapData) => {
                return new Promise((resolve) => {
                  /* Datapackage (-1 delemiter)
                            [0]	 ->	 fgsData [fgsTypeID, x * (LigandID, ModelID), -1, ...]
                        */
                  let check = FGSMapData["fgsMapData_recv"]["0"];
                  if (state.fgsMapData_recv === null && check != null) {
                    //console.log("FGSmapData", FGSMapData["fgsMapData_recv"]["0"]);
                    commit("setFGSMapData_recv", FGSMapData);
                    commit(
                      "setCurSelected_FGS_clusterData",
                      parseFGSMapData(check)
                    );
                    commit("setNew_FGSclusterData", true);
                  } else if (check != null) {
                    if (
                      state.fgsMapData_recv["fgsMapData_recv"]["0"].length !==
                        check.length ||
                      state.fgsMapData_recv["fgsMapData_recv"]["0"][
                        check.length - 1
                      ] !== check[check.length - 1]
                    ) {
                      //console.log("state.curSelected_FGS_clusterData.length !== data.length",state.curSelected_FGS_clusterData.length ," ", data.length);
                      // console.log("FGSmapData", FGSMapData["fgsMapData_recv"]["0"]);
                      commit("setFGSMapData_recv", FGSMapData);
                      commit(
                        "setCurSelected_FGS_clusterData",
                        parseFGSMapData(check)
                      );
                      commit("setNew_FGSclusterData", true);
                    }
                  }
                  clearInterval(interval);
                  setInterval(f1(), state.refreshTime);
                  console.log(
                    "refreshTime",
                    state.refreshTime,
                    "testMode",
                    state.testMode
                  );
                  resolve();
                });
              });
            //.then(() => console.log("END: parsing FGSMapData"));
          },

          state.refreshTime
        );
      }
      f1();
    },

    pollGUIData({ state, commit }) {
      function f1() {
        let interval = setInterval(
          () => {
            return fetch("/getGUIData", {
              method: "POST",
              headers: {
                "Content-Type": "application/json",
              },
            })
              .then((response) => response.json())
              .then((GUIdat) => {
                return new Promise((resolve) => {
                  /* Datapackage
                          [0]	 ->	 paramName;paramName....
                          [1]	 ->	 [floatVal, floatVal,...]
                        */

                  if (GUIdat["GUIdata_recv"]["1"][0] != state.GUIdataHash) {
                    commit("setBlockUserInput", false);
                    commit("setGUIdataHash", GUIdat["GUIdata_recv"]["1"][0]);
                    let GUIall = {};
                    let GUI_boolNames = [];
                    let GUI_floats = {};
                    let GUI_ints = {};
                    let GUIparamNames = GUIdat["GUIdata_recv"]["0"].split(";");
                    for (let i = 0; i < GUIparamNames.length; i++) {
                      let curParamName = GUIparamNames[i];
                      let curParamVal = GUIdat["GUIdata_recv"]["1"][i];
                      if (
                        curParamVal == 1.0 &&
                        !state.GUIIntParams.includes(curParamName) &&
                        !state.GUIfloatParams.includes(curParamName)
                      ) {
                        GUI_boolNames.push(curParamName);
                      }
                      if (state.GUIfloatParams.includes(curParamName)) {
                        GUI_floats[curParamName] = parseFloat(curParamVal);
                      }
                      if (state.GUIIntParams.includes(curParamName)) {
                        GUI_ints[curParamName] = parseInt(curParamVal);
                      }
                      //GUI[GUIparamNames[i]] = GUIdat["GUIdata_recv"]["1"][i];
                    }
                    GUIall["boolNames"] = GUI_boolNames;
                    GUIall["float"] = GUI_floats;
                    GUIall["int"] = GUI_ints;
                    //console.log("GUIdat", GUIdat);
                    //console.log("GUI:", GUIall);
                    commit("setGUIdata", GUIall);
                  }
                  clearInterval(interval);
                  setInterval(f1(), state.refreshTime);
                  resolve();
                });
              });
            //.then(() => console.log("END: parsing FGSMapData"));
          },

          state.refreshTime
        );
      }
      f1();
    },

    pollSelectionData({ state, commit }) {
      function f1() {
        let interval = setInterval(
          () => {
            return fetch("/getSelectionData", {
              method: "POST",
              headers: {
                "Content-Type": "application/json",
              },
            })
              .then((response) => response.json())
              .then((SelectDat) => {
                return new Promise((resolve) => {
                  /* Datapackage
                          [0]	 ->	 paramName;paramName....
                          [1]	 ->	 [floatVal, floatVal,...]
                        */
                  //console.log("SelectDat", SelectDat);
                  let MegaMol_selectionDataHash =
                    SelectDat["selectionData_recv"][1][0];
                  let clusterID = SelectDat["selectionData_recv"][1][1];
                  let ligandID = SelectDat["selectionData_recv"][1][2];
                  let modelID = SelectDat["selectionData_recv"][1][3];
                  //console.log("clusterID", clusterID)
                  //console.log("ligandID", ligandID)
                  //console.log("modelID", modelID)
                  if (
                    parseInt(MegaMol_selectionDataHash) !=
                    parseInt(state.MegaMol_selectionDataHash)
                  ) {
                    commit("setBlockUserInput", false);
                    commit("setSupressSendData", true);
                    if (clusterID == -1) {
                      if (store.state.showNoise == true) {
                        commit("setClusterID", "All");
                      } else {
                        commit("setClusterID", "All_withoutNoise");
                      }
                    } else {
                      commit("setClusterID", clusterID);
                    }
                    commit("setLigandID", ligandID);
                    commit("setModelID", modelID);
                    commit(
                      "setMegaMol_selectionDataHash",
                      MegaMol_selectionDataHash
                    );
                  }
                  clearInterval(interval);
                  setInterval(f1(), state.refreshTime);
                  resolve();
                }).then(() => {
                  commit("setSupressSendData", false);
                });
              });
            //.then(() => console.log("END: parsing FGSMapData"));
          },

          state.refreshTime
        );
      }
      f1();
    },
  },
  mutations: {
    setInfoOverlayContent(state, val) {
      let header = val[0];
      let filename = val[1];
      fetch("/getMarkdownFile/" + filename, {
        method: "GET",
        headers: {
          "Content-Type": "application/json",
        },
        //body: JSON.stringify({exportData: this.slider1.val})
      }).then((response) =>
        response.json().then((content) => {
          state.infoOverlayContent = [header, marked(content)];
        })
      );
    },
    getProLigInterPNG(state, val) {
      if (state.curProLigInterPNG_name == val[0]) {
        state.curProLigInterPNG = "";
        state.curProLigInterPNG_name = "";
      } else {
        state.curProLigInterPNG_name = val[0];

        console.log("ProLig: ", val);
        let existingKeys = Object.keys(state.proLigInterPNGs);
        let nonExistingKeys = val.filter((x) => !existingKeys.includes(x));
        if (nonExistingKeys.length > 0) {
          fetch("/getProLigInterPNG", {
            method: "POST",
            headers: {
              "Content-Type": "application/json",
            },
            body: JSON.stringify(nonExistingKeys),
          })
            .then((response) => response.json())
            .then(
              (newPNGs) =>
                (state.proLigInterPNGs = {
                  ...state.proLigInterPNGs,
                  ...newPNGs,
                })
            )
            .then(() => {
              state.curProLigInterPNG = state.proLigInterPNGs[val[0]];
            });
        } else {
          state.curProLigInterPNG = state.proLigInterPNGs[val[0]];
        }
      }
    },
    addSvgs(state, val) {
      let existingKeys = Object.keys(state.zincSVGs);
      let nonExistingKeys = val.filter((x) => !existingKeys.includes(x));
      fetch("/getZINCsvgs", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(nonExistingKeys),
      })
        .then((response) => response.json())
        .then(
          (newSvgs) => (state.zincSVGs = { ...state.zincSVGs, ...newSvgs })
        );
    },
    add_proLigInter_PNG(state, val) {
      fetch("/getProLigInterPNG", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(val),
      })
        .then((response) => response.json())
        .then(
          (newPNGs) =>
            (state.proLigInterPNGs = { ...state.proLigInterPNGs, ...newPNGs })
        );
    },
    setNrgDict(state, data) {
      state.nrgDict = data;
    },
    setInfoOverlayState(state, val) {
      state.showInfoOverlay = val;
    },
    setBlockUserInput(state, val) {
      console.log(
        "state.showBlockUserInput",
        val,
        "testMode: ",
        state.testMode
      );
      if (!state.testMode) {
        state.showBlockUserInput = val;
      }
      if (val) {
        console.log("wait to unlock...");
        delay(5000).then(() => {
          if (val) {
            state.showBlockUserInput = false;
          }
        });
      }
    },
    setStackingDat(state, data) {
      state.stackingDat = data;
    },
    setShowNavDrawer(state, val) {
      state.showNavDrawer = val;
    },
    setStackedBar_coloring(state, val) {
      state.stackedBar_coloring = val;
    },
    setCheckedFGSgroups(state, data) {
      state.checkedFGSgroups = data;
    },
    setFGSMapData_recv(state, data) {
      state.fgsMapData_recv = data;
    },
    setGUIcolorScaleSVGs(state, data) {
      state.GUIcolorScaleSVGs[data[0]] = data[1];
    },
    setCurSelected_FGS_clusterData(state, data) {
      state.curSelected_FGS_clusterData = data;
    },
    setNew_FGSclusterData(state, val) {
      state.new_FGSclusterData = val;
    },
    setClusterDataPrepared(state, data) {
      state.clusterData_prepared = data;
    },
    setOperationMode(state, data) {
      state.testMode = data;
    },
    setPollTime(state, data) {
      state.refreshTime = data;
    },
    setSelectedFuncGroup(state, val) {
      state.selectedFuncGroup = val == state.selectedFuncGroup ? null : val;
      console.log("val", val);
    },
    setColorMapSVGs(state, val) {
      state.colorMapSVGs = val;
    },
    setLgdMdl(state, val) {
      state.LgdMdl = val;
    },
    setGUIdataHash(state, val) {
      state.GUIdataHash = val;
    },
    setGUIdata(state, val) {
      state.GUIdata = val;
    },
    setMegaMol_selectionDataHash(state, val) {
      state.MegaMol_selectionDataHash = val;
    },
    setSupressSendData(state, val) {
      if (val == false) {
        setTimeout(function () {
          //console.log("wait to set suppress");
        }, 1500);
      }
      state.supressSendData = val;
    },
    setCurrLigandFuncDat(state, ligandID) {
      //console.log("ligandIDTTTT",ligandID);
      state.currLigandFuncDat =
        state.funcGroupData_recv[
          state.clusterData_prepared.lgdNameID[ligandID]
        ];
      // add PNG
      let existingKeys = Object.keys(state.functionalGroupsPNGs);
      let nonExistingKeys = Object.keys(store.state.currLigandFuncDat).filter(
        (x) => !existingKeys.includes(x)
      );
      fetch("/getLigandPNGs", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(nonExistingKeys),
      })
        .then((response) => response.json())
        .then(
          (newPngs) =>
            (state.functionalGroupsPNGs = {
              ...state.functionalGroupsPNGs,
              ...newPngs,
            })
        );
    },
    setClusterData(state, val) {
      state.clusterData_recv = val;
      /*
      state.maxEnergy = Math.ceil(d3.max(val[clusterDataKeys.energies].map(function (d){return Math.abs(d);})));
      state.minEnergy = Math.floor(d3.min(val[clusterDataKeys.energies].map(function (d, i){
        if (state.showNoise) {
          return Math.abs(d);
        } else {
				if (val[clusterDataKeys.clusterAssignment][i] == -1) {
					return state.max;
				} else {
					return Math.abs(d);
				}
			}})));
			*/
    },
    setSorting(state, val) {
      state.sorting = val;
    },
    setDynamicEnergyMin(state, val) {
      state.dynamicEnergyMin = val;
    },
    setFuncGroupData(state, data) {
      /* Datapackage
        [0]	->	funcGroupsPerLigands
				[LigandID, funcGroupsCnt, X * (funcGroupID, groupCnt), -1, LigandID, ...]
				-1 works as additional delimiter between ligands
        [1]	->	funcGroupsPerClusters
				[ClusterID, funcGroupsCnt, X * (funcGroupID, groupCnt, clusterCnt, Y * clusterSizes, -1), -2, ClusterID, ...]

		*** additional information ***
        funcGroupsCnt:	count of found functional groups
        X:				count of found functional groups
        funcGroupID:	specific ID for a type of functional group
        groupCnt:		count of how many times a specific functional group was found
        clusterCnt:		number of differentiable positions in a pocket where a specific func. groups is located
        Y:				number of differentiable positions in a pocket where a specific func. groups is located
        clusterSizes:	count of a specific func. group on such a position (see clusterCnt)
		-1:				is used as delimiter between ligands or func. groups
		-2:				is used as delimiter between clusters
    */

      let keys = funcGroupDataKeys;
      let perCluster = data[keys.funcGroupsPerClusters];
      let perLigand = data[keys.funcGroupsPerLigands];

      /******************************************
       ******* parse funcGroupsPerLigands *******
       ******************************************/

      //console.log("perCluster",perCluster)
      //console.log("perLigand",perLigand)

      let startIdx = 0;
      let idx = perLigand.indexOf(-1);
      let ligandSlices = [];
      let allLigFuncData = {};
      let invertedLigFuncData = {};

      while (idx != -1) {
        ligandSlices.push(perLigand.slice(startIdx, idx));
        startIdx = idx + 1;
        idx = perLigand.indexOf(-1, idx + 1);
      }

      ligandSlices.forEach((ligandArray) => {
        let ligandData = {};

        let zincName =
          state.clusterData_recv[clusterDataKeys.zincNames][ligandArray[0]];
        // build ligandData
        for (let i = 0; i < ligandArray[1]; i++) {
          ligandData[ligandArray[2 + i * 2]] = ligandArray[2 + i * 2 + 1];
          if (!(ligandArray[2 + i * 2] in invertedLigFuncData)) {
            invertedLigFuncData[ligandArray[2 + i * 2]] = [zincName];
          } else {
            invertedLigFuncData[ligandArray[2 + i * 2]].push(zincName);
          }
        }
        allLigFuncData[zincName] = ligandData;
      });

      state.funcGroupData_recv = allLigFuncData;
      state.funcGroupData_inverted = invertedLigFuncData;
      //console.log("INVERTED_FUNC_DICT ",Object.keys(state.funcGroupData_inverted).length, state.funcGroupData_inverted)

      /*****************************************
       ****** parse funcGroupsPerClusters ******
       *****************************************/
      /*
          [1]	->	funcGroupsPerClusters
          [ClusterID, funcGroupsCnt, X * (funcGroupID, groupCnt, groupCnt * [liganID,modelID], -1), -2, ClusterID, ...]
      */
      let allClusterFuncData = [];
      let dummyData_perFuncGroupData = [];
      let dummyData_allPerClusterFuncData = [];
      startIdx = 0;
      idx = 0;

      // initialize
      dummyData_perFuncGroupData["funcGroupID"] = 0;
      dummyData_perFuncGroupData["funcGroupCnt"] = 0;
      dummyData_perFuncGroupData["lig_and_mdl_IDs"] = 0;
      for (let j = 0; j < state.numberOf_FGStypeIDs; j++) {
        dummyData_allPerClusterFuncData[j] = Object.assign(
          {},
          dummyData_perFuncGroupData
        );
      }

      // per Cluster Data
      while (idx < perCluster.length) {
        if (perCluster[idx] == -1 && idx != 0) {
          idx++;
        } else {
          if (idx != 0) {
            console.log("missing clusterData delimiter");
          }
        }
        if (idx >= perCluster.length) {
          break;
        }

        if (perCluster[idx] == -2 && idx != 0) {
          idx++;
        } else {
          if (idx != 0) {
            console.log("missing clusterData delimiter2");
          }
        }

        let clusterID = perCluster[idx + 0];
        let foundFuncGroupsCnt = perCluster[idx + 1];
        idx += 2;
        // per FuncGroup Data

        let allPerClusterFuncData = [];
        allPerClusterFuncData = Object.assign(
          {},
          dummyData_allPerClusterFuncData
        );

        let z = idx;
        for (let i = 0; i < foundFuncGroupsCnt; i++) {
          let perFuncGroupData = [];

          if (perCluster[z] == -1 && i != 0) {
            z++;
          } else {
            if (i != 0) {
              console.log("missing funcGroup delimiter");
            }
          }

          let funcGroupID = perCluster[z];
          let funcGroupCnt = perCluster[z + 1];

          perFuncGroupData["funcGroupID"] = funcGroupID;
          perFuncGroupData["funcGroupCnt"] = funcGroupCnt;

          let y;
          let ligMdlIDs = [];
          let pairsOf = 2;
          for (let j = 0; j < funcGroupCnt; j++) {
            y = z + 2 + j * pairsOf;
            // push ligandID
            ligMdlIDs[j * pairsOf + 0] = perCluster[y + 0];
            // push modelID
            ligMdlIDs[j * pairsOf + 1] = perCluster[y + 1];
          }

          perFuncGroupData["lig_and_mdl_IDs"] = ligMdlIDs;

          // push funcGroupData
          allPerClusterFuncData[funcGroupID] = perFuncGroupData;
          z = y + pairsOf;
        }
        idx = z;

        // push clusterData
        //console.log("allPerClusterFuncData",allPerClusterFuncData);
        allClusterFuncData[clusterID] = allPerClusterFuncData;
      }

      //console.log("allClusterFuncData",allClusterFuncData);
      state.clusterFuncData = allClusterFuncData;
    },
    setLargeLigand(state, val) {
      state.largeLigand = val;
    },
    setShowNoise(state, val) {
      state.showNoise = val;
      if (state.clusterID == "All") {
        state.clusterID = "All_withoutNoise";
      } else if (state.clusterID == "All_withoutNoise") {
        state.clusterID = "All";
      }
    },
    setCheckB_fgsClust(state, val) {
      state.checkB_fgsClust = val;
      console.log("checkB_fgsClust", state.checkB_fgsClust);
    },
    setTableHeaders(state, val) {
      state.tableHeaders = val;
    },
    setTableData(state, val) {
      state.tableData = val;
    },
    setTableSearch(state, val) {
      state.tableSearch = val;
    },
    setTableHeaders2(state, val) {
      state.tableHeaders2 = val;
    },
    setTableData2(state, val) {
      state.tableData2 = val;
    },
    setSubmitVal(state, val) {
      state.submitVal = val;
    },
    setLgdEnergies(state, val) {
      state.lgdEnergies = val;
    },
    setMaxEnergy(state, val) {
      state.maxEnergy = val;
    },
    setMinEnergy(state, val) {
      state.minEnergy = val;
    },
    calcMaxEnergy(state) {
      state.maxEnergy = Math.ceil(
        d3.max(
          state.clusterData_recv[clusterDataKeys.energies].map(function (d) {
            return Math.abs(d);
          })
        )
      );
    },
    calcMinEnergy(state) {
      state.minEnergy = Math.floor(
        d3.min(
          state.clusterData_recv[clusterDataKeys.energies].map(function (d, i) {
            if (Math.abs(d) != 0) {
              if (state.showNoise) {
                return Math.abs(d);
              } else {
                if (
                  state.clusterData_recv[clusterDataKeys.clusterAssignment][
                    i
                  ] == -1
                ) {
                  return state.maxEnergy;
                } else {
                  return Math.abs(d);
                }
              }
            }
          })
        )
      );
    },
    setGblMdlID(state, gblMdlID) {
      state.gblMdlID = parseInt(gblMdlID);
    },
    setLigandID(state, ligandID) {
      state.ligandID = parseInt(ligandID);
      state.selectionDataHash++;
    },
    setModelID(state, modelID) {
      state.modelID = parseInt(modelID);
      state.selectionDataHash++;
    },
    setLigandModelID(state, ligandModelID) {
      state.ligandID = parseInt(ligandModelID[0]);
      state.modelID = parseInt(ligandModelID[1]);
      state.selectionDataHash++;
    },
    setClusterID(state, clusterID) {
      state.clusterID = clusterID;
      state.selectionDataHash++;
    },
    setMdlInSelClust(state, val) {
      state.mdlsInSelectedClust = val;
    },
    setIsFullscreen(state, val) {
      console.log("isFullscreen", val);
      state.isFullscreen = val;
    },
    setLgdTableTitle(state, ligandID) {
      state.lgdTableTitle =
        "(" + ligandID + ") " + state.clusterData_prepared.lgdNameID[ligandID];
      if (
        typeof state.clusterData_prepared.lgdNameID[ligandID] === "undefined" ||
        ligandID === -1
      ) {
        state.lgdTableTitle = "no selection";
      }
    },
    setPresentClustData(state, val) {
      state.presentClustData = val;
    },
    setHbonds(state, val) {
      state.hbonds = val;
    },
    setHalogenBonds(state, val) {
      state.halogenBonds = val;
    },
    setHydrophobicInteractions(state, val) {
      state.hydrophobicInteractions = val;
    },
    setMetalComplexes(state, val) {
      state.metalComplexes = val;
    },
    setPiCationInteractions(state, val) {
      state.piCationInteractions = val;
    },
    setPiStacks(state, val) {
      state.piStacks = val;
    },
    setSaltBridges(state, val) {
      state.saltBridges = val;
    },
    setEfficiency(state, val) {
      state.efficiency = val;
    },
    setHighlightBar(state, val) {
      state.highlightBar = val;
    },
    setBarScaleX(state, val) {
      state.barScaleX = val;
    },
    setBarScaleY(state, val) {
      state.barScaleY = val;
    },
    setTreeActives(state, val) {
      state.treeActives = val;
    },
  },
});
new Vue({
  vuetify,
  store: store,
  render: (h) => h(App),
}).$mount("#app");
