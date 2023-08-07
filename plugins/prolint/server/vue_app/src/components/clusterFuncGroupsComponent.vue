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

    <v-card-title class="pt-4">
      Functional Groups:
      {{
        clusterID == "All"
          ? "All Clusters with noise "
          : clusterID == "All_withoutNoise"
          ? "All Clusters"
          : "Cluster " + clusterID
      }}
      <infoOverlayComponent header="Functional Groups View:">
      </infoOverlayComponent>
      <v-spacer></v-spacer>
      <v-row>
        <v-col cols="4"></v-col>
        <v-col>
          <v-select
            class="pb-0"
            :items="filterModeOptions"
            v-model="checkB_fgsClust"
            label="filter mode"
            dense
            hint=" filter mode - Cluster Statistics"
            persistent-hint
            single-line
            flat
          ></v-select>
        </v-col>
      </v-row>
    </v-card-title>
    <v-card id="clusterTreeView">
      <v-treeview
        dense
        hoverable
        color="#0e8cf3"
        multiple-active
        :selectable="SetSelectable"
        v-model="selectionData"
        :items="treeviewItems"
        :active="setActives"
        return-object
      >
        <template v-slot:append="{ item }">
          <div>
            {{
              item.part != null ? item.part + "/" + item.amount : item.amount
            }}
          </div>
        </template>
      </v-treeview>
    </v-card>

    <!--
    <div id ="funcGroupTooltip">  </div>
    <div id ="functionalGroupsHorizonalAxis"></div>
    <div id ="scrollWrappper">
      <div id ="functionalGroupsMain"></div>
    </div>
    -->
  </v-container>
</template>

<script>
//import {functionalGroupChart} from "@/assets/js/functionalGroupChart"
import { mapState } from "vuex";
import fgClassifications from "@/assets/js/fgClassifications";
import { clusterDataKeys as k } from "../main";
import { updateLgdTableData } from "@/assets/js/cluster";
import infoOverlayComponent from "./infoOverlayComponent";
import { scrollParentToChild } from "./componentUtils";

export default {
  // name of the component
  name: "clusterFuncGroupsComponent",
  components: {
    infoOverlayComponent,
  },
  // data section of the Vue component. Access via this.<varName> .
  data: () => ({
    selectionData: [],
    clusterFuncCounts: {},
    loading: false,
    sorting: null,
    show_noise: false,
    funcDataChart: null,
    firstChart: -1,
  }),
  computed: {
    ...mapState({
      clusterID: (state) => state.clusterID,
      clusterFuncData: (state) => state.clusterFuncData,
      functionalGroupsIDsToWord: (state) => state.functionalGroupsIDsToWord,
      funcGroupData_inverted: (state) => state.funcGroupData_inverted,
      clusterData_recv: (state) => state.clusterData_recv,
      checkB_fgsClust: (state) => state.checkB_fgsClust,
      curSelected_FGS_clusterData: (state) => state.curSelected_FGS_clusterData,
      treeActives: (state) => state.treeActives,
      new_FGSclusterData: (state) => state.new_FGSclusterData,
      scrollToRow: (state) => state.curSelected_FGS_clusterData,
    }),
    filterModeOptions: {
      get: function () {
        return [
          { text: "Functional Group Cluster", value: true },
          { text: "Manual Treeview Selection", value: false },
        ];
      },
      set: function () {
        if (this.new_FGSclusterData) {
          return true;
        }
      },
    },
    SetSelectable: {
      get: function () {
        if (this.checkB_fgsClust) {
          return false;
        } else {
          return true;
        }
      },
      set: function () {},
    },
    checkB_fgsClust: {
      get() {
        return this.$store.state.checkB_fgsClust;
      },

      set(val) {
        //this.$store.commit("setCheckB_fgsClust", true);
        if (val == true) {
          this.$store.commit("setNew_FGSclusterData", true);
          let GUI = this.$store.state.GUIdata;
          let val = ["funcGroupCentroids", 1];
          GUI["boolNames"][val[0]] = val[1];
          this.$store.commit("setGUIdata", GUI);
          this.$store.commit("setBlockUserInput", true);
          fetch("/navigation", {
            method: "POST",
            headers: {
              "Content-Type": "application/json",
            },
            body: JSON.stringify({ exportData: val }),
          });
        } else {
          //console.log("this.setSelection(null)");
          this.setSelection([]);
        }
        this.$store.commit("setCheckB_fgsClust", val);
      },
    },

    setActives: {
      get: function () {
        if (this.checkB_fgsClust == true && this.new_FGSclusterData == true) {
          //console.log("here in set active...");
          //this.active = this.treeActives;
          this.setSelection(this.treeActives);
          if (this.selectionData.length == this.treeActives.length) {
            // TODO: make this work
            //this.$store.commit("setNew_FGSclusterData",false);
          }

          return this.treeActives;
        } else {
          this.setSelection(this.selectionData);
          return this.selectionData;
        }
      },
    },

    treeviewItems: {
      //TODO: Specificaiton tells us that objects are unordered atm this func relies on order
      get: function () {
        if (this.clusterFuncData != null) {
          if (this.clusterID == null) {
            this.$store.commit("setClusterID", "All");
          }

          /*****************************************
           **** TREE VIEW - PARSE FUNC GROUPS ******
           *****************************************/
          this.$store.commit("setTreeActives", []);
          //console.log("CLUSTER TreeviewPart")
          let outArr = [];
          let id = 1;
          // eslint-disable-next-line vue/no-side-effects-in-computed-properties
          this.clusterFuncCounts = {};
          if (this.clusterID == "All" || this.clusterID == "All_withoutNoise") {
            for (const [clstrID] of Object.entries(this.clusterFuncData)) {
              for (const [element, value] of Object.entries(
                this.clusterFuncData[clstrID]
              )) {
                if (value["funcGroupCnt"] > 0) {
                  if (element in this.clusterFuncCounts) {
                    // eslint-disable-next-line vue/no-side-effects-in-computed-properties
                    this.clusterFuncCounts[element] += value["funcGroupCnt"];
                  } else {
                    // eslint-disable-next-line vue/no-side-effects-in-computed-properties
                    this.clusterFuncCounts[element] = value["funcGroupCnt"];
                  }
                }
              }
            }
          } else {
            // dict. of funcGroupIDs their Counts:  [funcGroupID] = Count
            for (const [element, value] of Object.entries(
              this.clusterFuncData[this.clusterID]
            )) {
              // eslint-disable-next-line vue/no-side-effects-in-computed-properties
              this.clusterFuncCounts[element] = value["funcGroupCnt"];
            }
          }
          //console.log("CLUSTER FuncGroups: ",this.clusterFuncCounts)
          let clusterFuncCountsCopy = Object.assign({}, this.clusterFuncCounts);

          // invert order of funcGroupIDs-Count-dictionary: now big to little
          let clusterFuncCountsKeys = Object.keys(
            clusterFuncCountsCopy
          ).reverse();

          // execute until all Counts in funcGroupIDs-Count-dictionary are 0
          while (Object.values(clusterFuncCountsCopy).some((e) => e > 0)) {
            clusterFuncCountsKeys.forEach((element) => {
              while (clusterFuncCountsCopy[element] > 0) {
                let tarArray = [];
                this.functionalGroupRecursion(tarArray, element);
                //console.log(tarArray)
                tarArray.forEach((elementInner) => {
                  clusterFuncCountsCopy[elementInner]--;
                });
                outArr.push(this.formatTreeviewEntry(tarArray, id));
                id += tarArray.length;
              }
            });
          }
          //for outermost part
          outArr = this.mergeDuplicates(outArr);
          //for rest
          outArr = this.recursiveMerge(outArr);
          //console.log(outArr)
          outArr = this.addUndefinedGroups(outArr);
          // set the id list of the treeView nodes for that fgsTypes which are contained in fgsClustData
          outArr.forEach((element) => {
            //this.treeActives.push(element);
            if (
              this.checkB_fgsClust &&
              this.curSelected_FGS_clusterData != null
            ) {
              let tmp = Object.assign(
                {},
                this.curSelected_FGS_clusterData["fgsTypeIDs"]
              );
              if (
                Object.values(tmp).indexOf(parseInt(element.fgsTypeID)) > -1
              ) {
                this.treeActives.push(element);
              }
            }
          });
          //console.log("TreeView",outArr);
          //console.log("this.treeActives",this.treeActives);
          return outArr;
        }
        return [];
      },
    },
  },

  watch: {
    scrollToRow: function () {
      //https://www.appsloveworld.com/vuejs/100/76/scroll-to-a-programmatically-selected-row-in-a-v-data-table
      /** @noinspection */
      setTimeout(() => {
        const row = document
          .getElementById("clusterTreeView")
          .getElementsByClassName("v-treeview-node--selected")[0];

        if (row) {
          const parent = document.getElementById("clusterTreeView");
          scrollParentToChild(parent, row);
        }
      }, 10);
    },
    selectionData: function () {
      let filterZINC = [];
      this.$store.commit("setCheckedFGSgroups", this.selectionData);
      updateLgdTableData(0);
      this.selectionData.forEach((element) => {
        if (this.clusterID == "All" || this.clusterID == "All_withoutNoise") {
          for (const [clstrID] of Object.entries(this.clusterFuncData)) {
            let fgsCount = 0;
            let iter = 0;
            let value = this.clusterFuncData[clstrID][element.fgsTypeID];
            if (value["funcGroupCnt"] > 0) {
              fgsCount = value["funcGroupCnt"];
              iter = value["lig_and_mdl_IDs"];
              for (let i = 0; i < fgsCount; i++) {
                let ligID = iter[i * 2];
                if (this.checkB_fgsClust) {
                  // only push ligID if contained in func. group cluster data
                  //console.log("ligID",ligID)
                  if (
                    this.curSelected_FGS_clusterData[element.fgsTypeID] != null
                  ) {
                    if (
                      Object.values(
                        this.curSelected_FGS_clusterData[element.fgsTypeID][
                          "ligID"
                        ]
                      ).indexOf(parseInt(ligID)) > -1
                    ) {
                      filterZINC.push(
                        this.clusterData_recv[k.zincNames][ligID]
                      );
                    }
                  }
                } else {
                  filterZINC.push(this.clusterData_recv[k.zincNames][ligID]);
                }
              }
            }
          }
        } else {
          //console.log("element.fgsTypeID", element.fgsTypeID);
          //filterZINC.push(...this.funcGroupData_inverted[element.fgsTypeID]);
          let fgsCount =
            this.clusterFuncData[this.clusterID][element.fgsTypeID][
              "funcGroupCnt"
            ];
          let iter =
            this.clusterFuncData[this.clusterID][element.fgsTypeID][
              "lig_and_mdl_IDs"
            ];
          for (let i = 0; i < fgsCount; i++) {
            let ligID = iter[i * 2];
            if (this.checkB_fgsClust) {
              // only push ligID if contained in func. group cluster data
              //console.log("ligID", ligID);
              if (this.curSelected_FGS_clusterData[element.fgsTypeID] != null) {
                if (
                  Object.values(
                    this.curSelected_FGS_clusterData[element.fgsTypeID]["ligID"]
                  ).indexOf(parseInt(ligID)) > -1
                ) {
                  filterZINC.push(this.clusterData_recv[k.zincNames][ligID]);
                }
              }
            } else {
              filterZINC.push(this.clusterData_recv[k.zincNames][ligID]);
            }
          }
        }
      });
      let uniqueItems = [...new Set(filterZINC)];
      this.$store.commit("setTableSearch", uniqueItems.join(","));
    },
    /*
    clusterID: function() {
      console.log("CLUSTERID", this.clusterID)
      if (this.firstChart == -1){
        console.log(this.clusterFuncData[this.clusterID])
        this.funcDataChart = new functionalGroupChart(["#functionalGroupsHorizonalAxis","#functionalGroupsMain"], this.clusterFuncData[this.clusterID])
        this.firstChart = 0
      }
      else{
        this.funcDataChart.updateData(this.clusterFuncData[this.clusterID])
      }

    }
    */
  },

  // functions to call on mount (after DOM etc. is built)
  mounted() {
    //this.$vuetify.theme.themes.light.primary = '#2196F3';
  },

  // methods of this vue component access via this.<funcName>
  methods: {
    formatTreeviewEntry(unformattedArray, id) {
      let outObj = 0;
      for (let index = 0; index < unformattedArray.length; index++) {
        const element = unformattedArray[index];
        let part = null;
        if (this.checkB_fgsClust && this.curSelected_FGS_clusterData != null) {
          let tmp = Object.assign(
            {},
            this.curSelected_FGS_clusterData["fgsTypeIDs"]
          );
          if (Object.values(tmp).indexOf(parseInt(element)) > -1) {
            part = this.curSelected_FGS_clusterData[element]["ligID"].length;
          }
        }

        if (outObj == 0) {
          outObj = {
            id: id,
            fgsTypeID: element,
            name: this.functionalGroupsIDsToWord[element],
            amount: parseInt(this.clusterFuncCounts[element]),
            part: part,
            children: [],
          };
        } else {
          let tempObj = {
            id: id + unformattedArray.length - index,
            fgsTypeID: element,
            amount: parseInt(this.clusterFuncCounts[element]),
            part: part,
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
            if (element.fgsTypeID == compareElem.fgsTypeID) {
              if (element.children) {
                element.children = element.children.concat(
                  compareElem.children
                );
                /*
                let amtSubgroups = this.sumAmountSubgroups(element.children)
                //console.log("amtSubgrp", amtSubgroups)
                let diffAmount = element.amount - amtSubgroups
                //console.log("diffAmount", diffAmount)

                if(diffAmount >0){
                  element.children.push({id:-1, name: "not further specified", amount: diffAmount, children: []})
                }
                */
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

    sumAmountSubgroups(subgroups) {
      let outAmount = 0;
      subgroups.forEach((element) => {
        //console.log("elemAmt", element.amount)
        outAmount += element.amount;
      });
      return outAmount;
    },

    addUndefinedGroups(array) {
      array.forEach((element) => {
        if (element.children) {
          //console.log("RECURSIVE MERGE STEP: ", element)
          if (element.children.length > 0) {
            this.addUndefinedGroups(element.children);
          }
          let pushElem = this.undefinedGroupsRecursion(element);

          if (Object.keys(pushElem).length > 0 && element.children.length > 0) {
            //console.log("pushElem",pushElem);
            element.children.push(pushElem);
          }
        }
      }, array);
      // console.log("undef added", array);
      return array;
    },

    undefinedGroupsRecursion(elem) {
      let amtSubgroups = this.sumAmountSubgroups(elem.children);
      //console.log("amtSubgrp", amtSubgroups)
      let diffAmount = elem.amount - amtSubgroups;
      //console.log("diffAmount", diffAmount)

      if (diffAmount > 0) {
        return {
          id: elem.fgsTypeID + "_unspec",
          fgsTypeID: elem.fgsTypeID,
          name: "not further specified",
          amount: diffAmount,
          children: [],
        };
      } else {
        //console.log("no unspec.")
        return {};
      }
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
    setSelection(data) {
      // console.log("setSelection",data.length);
      this.selectionData = data;
    },
  },
};
</script>

<style>
@import "../assets/css/functionalGroupChart.css";
</style>
