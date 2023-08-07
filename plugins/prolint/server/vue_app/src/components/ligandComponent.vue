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

    <v-row>
      <!--- Larger Div --->
      <v-col cols="12" class="d-flex flex-column">
        <!-- was col = 8-->
        <v-card class="mx-auto my-auto colCard">
          <div id="LargeDiv">
            <div id="tooltipLigand"></div>
            <div id="ligand"></div>
            <v-row>
              <v-col class="pa-3" cols="5"></v-col>
              <v-col class="d-flex align-center pa-3" cols="2">
                <v-text-field
                  hide-details
                  label="~Bins (step 5)"
                  v-model="spinnerHisto"
                  class="mt-0 pt-0 pa-0"
                  type="number"
                  min="5"
                  max="100"
                  step="5"
                  outlined
                  style="width: 60px"
                ></v-text-field>
              </v-col>
            </v-row>
            <v-row></v-row>
          </div>
          <!--end LargeDiv-->
        </v-card>
      </v-col>
    </v-row>

    <v-row>
      <v-col cols="12" class="d-flex flex-column">
        <!-- was col = 8-->
        <v-card class="mx-auto my-auto colCard">
          <div id="scroll-wrapper-lig">
            <div class="card" id="smallMult"></div>
            <v-row v-for="(val, idx) in rowList" :key="idx" no-gutters>
              <v-col
                v-for="(val2, idx2) in val"
                :key="idx2"
                cols="4"
                :id="val2"
              >
              </v-col>
            </v-row>
          </div>
        </v-card>
      </v-col>
    </v-row>
  </v-container>
</template>

<script>
import {
  drawLigandHisto,
  plot_order,
  smallMultiples,
} from "@/assets/js/ligand";
import { mapState } from "vuex";

export default {
  name: "ligandComponent",
  data: () => ({
    loading: false,
    sampleTextInput: "SAMPLE",
    spinnerHisto: 10,
    ligandHistoData: null,
    plotOrder: [],
    rowList: [],
    smallMultiplesComponent: "",
  }),
  computed: {
    ...mapState({
      submitVal: (state) => state.submitVal,
      clusterData_recv: (state) => state.clusterData_recv,
      lgdMdl: (state) => state.lgdMdl,
    }),
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
    spinnerHisto: function () {
      document.getElementById("tooltipLigand").innerHTML = "";
      document.getElementById("ligand").innerHTML = "";
      document.getElementById("smallMult").innerHTML = "";
      this.ligandHistoData = drawLigandHisto(
        this.clusterData_recv,
        this.spinnerHisto,
        this.lgdMdl
      );
    },
    submitVal: function () {
      this.submit();
    },
    rowList: function () {
      this.$nextTick(() => {
        smallMultiples(this.plotOrder, this.ligandHistoData);
      });
    },
    clusterData_recv: function () {
      if (this.clusterData_recv != null) {
        this.ligandHistoData = drawLigandHisto(
          this.clusterData_recv,
          10,
          this.lgdMdl
        );
        this.plotOrder = plot_order(this.ligandHistoData.lgdEnergies);
        this.rowList = this.getColIds(this.plotOrder);
      }
    },
  },
  mounted() {
    if (!(this.clusterData_recv == null)) {
      this.ligandHistoData = drawLigandHisto(
        this.clusterData_recv,
        10,
        this.lgdMdl
      );
      this.plotOrder = plot_order(this.ligandHistoData.lgdEnergies);
      this.rowList = this.getColIds(this.plotOrder);
    }
  },
  // methods of this vue component access via this.<funcName>
  methods: {
    submit() {
      this.loading = true;
      fetch("/ligand", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({ exportData: this.submitVal }),
      })
        .then((response) => response.json())
        .then((dataContent) => {
          this.importedData = dataContent;
          console.log("Imported: ", this.importedData);
        })
        .finally(() => {
          this.loading = false;
        });
    },
    getColIds(plotOrder) {
      let outArray = [];
      let rowArray = [];
      for (let index = 0; index < plotOrder.length; index++) {
        if (index % 3 == 0 && index != 0) {
          outArray.push(rowArray);
          rowArray = [];
        }
        rowArray.push(`smallMultiples_${index}`);
      }
      return outArray;
    },
  },
};
</script>
