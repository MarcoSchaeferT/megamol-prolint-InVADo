<style>
@import "../assets/css/tooltipComponent.css";
</style>

<template>
  <v-scroll-y-transition>
    <div class="tooltipComponent" v-if="awesome" :style="getStyle()">
      <div class="my-2" v-show="awesome">
        <img :src="currentlySelectedFunctionalGroup" class="imgTooltip" />

        <v-btn @click="closeTooltip" color="black" fab dark id="closeButton">
          <v-icon>mdi-close</v-icon>
        </v-btn>
      </div>
    </div>
  </v-scroll-y-transition>
</template>

<script>
import { mapGetters, mapState } from "vuex";

export default {
  name: "tooltipSkeletalFormulaComponent",
  show: false,
  clientWidth: 0,
  clientHeight: 0,

  data: () => ({
    loading: false,
  }),
  computed: {
    ...mapState({
      selectedFuncGroup: (state) => state.selectedFuncGroup,
    }),
    ...mapGetters(["getFunctionalGroupPng"]),
    awesome: {
      get: function () {
        if (this.selectedFuncGroup != null) {
          return true;
        } else {
          return false;
        }
      },
    },
    currentlySelectedFunctionalGroup: {
      get: function () {
        //console.log("IMG FUNC LAUNCHED",this.selectedFuncGroup)
        if (this.selectedFuncGroup != null) {
          //console.log("COMPUTED TRIGGERED", this.selectedFuncGroup)
          let base64img = this.getFunctionalGroupPng(this.selectedFuncGroup);
          //console.log("image recieved", base64img)
          //let convertedImg = new Image()
          //convertedImg.src= base64img
          return base64img;
        }
        return null;
      },
    },
  },
  // Setup a Listener for one or more variables of "data", i.e. slider1 references to data -> slider1
  watch: {},

  // methods of this vue component access via this.<funcName>
  methods: {
    closeTooltip() {
      this.$store.commit("setSelectedFuncGroup", null);
    },
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
    getStyle() {
      let halfBigComponent =
        document.getElementById("clusterViewID")["clientHeight"];
      let clusterViewWidth =
        document.getElementById("clusterViewID")["clientWidth"];
      let halfBigComponent_all = document.getElementById(
        "dockingOverviewComponent_all"
      )["clientWidth"];

      let navDrawer = document.getElementById("navDrawerS");
      this.clientWidth = 360;
      clusterViewWidth -= this.clientWidth + 20;
      halfBigComponent -= this.clientWidth + 10;

      let navDrawerWidth = parseFloat(navDrawer["clientWidth"]); //+ halfBigComponent["clientWidth"] * 0.02;
      halfBigComponent_all += navDrawerWidth + clusterViewWidth;
      let css =
        "position: fixed; left: " +
        halfBigComponent_all +
        "px;  top: " +
        halfBigComponent +
        "px; width: " +
        this.clientWidth * 0.99 +
        "px; height: " +
        this.clientWidth * 0.99 +
        "px;";

      return css;
      //console.log("dockingOverviewComponent", tmp,width, height);
    },
  },
};
</script>
