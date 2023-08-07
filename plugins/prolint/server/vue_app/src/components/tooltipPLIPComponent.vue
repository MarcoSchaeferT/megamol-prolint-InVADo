<style>
@import "../assets/css/tooltipComponent.css";
</style>

<template>
  <v-scroll-y-transition>
    <div class="tooltipComponent" v-if="awesome" :style="getStyle()">
      <div class="my-2" v-show="awesome">
        <img :src="curProLigInterPNG" id="curLigPNG" />

        <v-btn
          @click="closeTooltip_PLIP"
          color="black"
          fab
          dark
          id="closeButton"
        >
          <v-icon>mdi-close</v-icon>
        </v-btn>
      </div>
    </div>
  </v-scroll-y-transition>
</template>

<script>
import { mapState } from "vuex";

export default {
  name: "tooltipPLIPComponent",
  show: false,
  clientWidth: 0,
  clientHeight: 0,

  data: () => ({
    loading: false,
  }),
  computed: {
    ...mapState({
      curProLigInterPNG: (state) => state.curProLigInterPNG,
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
  // Setup a Listener for one or more variables of "data", i.e. slider1 references to data -> slider1
  watch: {},

  // methods of this vue component access via this.<funcName>
  methods: {
    closeTooltip_PLIP() {
      this.$store.commit("getProLigInterPNG", []);
      console.log("closeTooltip_PLIP");
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
      let halfBigComponent = document.getElementById(
        "dockingOverviewComponent"
      );
      let halfBigComponent_all = document.getElementById(
        "dockingOverviewComponent_all"
      );
      let navDrawer = document.getElementById("navDrawerS");
      this.clientWidth = halfBigComponent["clientWidth"];
      let all = halfBigComponent_all["clientWidth"];
      let dif = all - this.clientWidth;
      let secDif = this.clientWidth * 0.001;
      let navDrawerWidth = parseFloat(navDrawer["clientWidth"]) + dif + secDif; //+ halfBigComponent["clientWidth"] * 0.02;
      let css =
        "position: fixed; left: " +
        navDrawerWidth +
        "px;  top: 10px; width: " +
        this.clientWidth * 0.99 +
        "px;";

      return css;
      //console.log("dockingOverviewComponent", tmp,width, height);
    },
  },
};
</script>
