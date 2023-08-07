<style>
@import "assets/css/main_css.css";
@import "assets/css/ligandStyles.css";
@import "assets/css/clusterStyles.css";
</style>
<template>
  <v-app>
    <v-overlay v-model="infoBoxSegement">
      <v-flex xs12 class="mx-auto my-auto pa-0 pb-0 colCard">
        <div class="infoComponent">
          <v-card color="blue-grey darken-2" class="white--text">
            <v-app-bar absolute scroll-target="infoOverlay">
              <v-card-title primary-title>
                <div class="headline">
                  {{ infoOverlayContentHeader }}
                </div>
              </v-card-title>
              <v-spacer></v-spacer>
              <v-btn @click="clickAction(false)" fab color="black" dark small>
                <v-icon>mdi-close</v-icon>
              </v-btn>
            </v-app-bar>
            <v-card-text>
              <div id="infoOverlay" class="overflow-y-auto">
                <div id="infoOverlay" v-html="infoOverlayContentText"></div>
              </div>
            </v-card-text>
            <!--
            <v-card-actions>
              <v-btn dark @click="clickAction(false)">close</v-btn>
            </v-card-actions>
            -->
          </v-card>
        </div>
      </v-flex>
    </v-overlay>

    <v-overlay v-model="showBlockUserInput">
      <v-flex xs12 class="mx-auto my-auto pa-0 pb-0 colCard">
        <div class="showBlockUserInputOverlay">
          <v-progress-circular
            :size="70"
            :width="7"
            color="gray"
            indeterminate
          ></v-progress-circular>
        </div>
      </v-flex>
    </v-overlay>
    <v-navigation-drawer
      app
      v-model="drawer"
      :expand-on-hover="isPinned"
      hide-overlay
      color="primary"
      dense
      width="512"
      id="navDrawerS"
      temporary
      fixed
      permanent
      absolute
    >
      <v-list-item class="px-2 pt-2">
        <v-btn elevation="2" fab small>
          <v-app-bar-nav-icon></v-app-bar-nav-icon>
        </v-btn>

        <div class="d-flex align-center">
          <v-avatar style="background-color: whitesmoke; left: 33px"
            ><v-spacer></v-spacer>
            <v-img
              class="shrink mr-2"
              contain
              :src="require('./assets/Pocket.svg')"
              transition="scale-transition"
              style="left: 3px; width: 100%; height: 140%"
            ></v-img>
          </v-avatar>

          <span
            style="
              position: absolute;
              left: 135px;
              width: 100%;
              font-size: x-large;
            "
            class="mr-2 white--text"
          >
            InVADo</span
          >
          <div style="position: fixed; right: 25px">
            <infoOverlayComponent
              class="pt-1"
              :icon-color="'white'"
              header="Sidebar Menu:"
            >
            </infoOverlayComponent>
          </div>
        </div>
      </v-list-item>
      <v-row>
        <v-col cols="1">
          <v-list-item class="px-2 pt-2">
            <v-tooltip right>
              <template v-slot:activator="{ on, attrs }">
                <v-btn
                  style="vertical-align: middle"
                  elevation="22"
                  fab
                  small
                  depressed
                  color="gray"
                  v-bind="attrs"
                  v-on="on"
                  @click.stop="isPinned = !isPinned"
                >
                  <v-icon color="black" v-if="!isPinned"
                    >mdi-pin-outline</v-icon
                  >
                  <v-icon color="grey darken-1" v-if="isPinned"
                    >mdi-pin-off-outline</v-icon
                  >
                </v-btn>
              </template>
              <span>auto-hide controls on/off</span>
            </v-tooltip>
          </v-list-item>
          <v-list-item class="px-2 pt-2 fillheight">
            <v-tooltip right>
              <template v-slot:activator="{ on, attrs }">
                <v-btn
                  style="vertical-align: middle"
                  elevation="22"
                  fab
                  small
                  depressed
                  color="gray"
                  v-bind="attrs"
                  v-on="on"
                  @click.stop="goFullscreen = !goFullscreen"
                >
                  <v-icon v-if="!isFullscreen" color="grey darken-1"
                    >mdi-fullscreen</v-icon
                  >
                  <v-icon v-if="isFullscreen" color="black"
                    >mdi-fullscreen-exit</v-icon
                  >
                </v-btn>
              </template>
              <span>fullscreen on/off</span>
            </v-tooltip>
          </v-list-item>
          <v-container fill-height><p></p></v-container>
        </v-col>
        <v-col>
          <v-list dense>
            <v-list-item>
              <v-list-item-content style="display: none"></v-list-item-content>
              <v-list-item-content>
                <navDrawer style="z-index: 0"></navDrawer>
              </v-list-item-content>
            </v-list-item>
          </v-list>
        </v-col>
      </v-row>
    </v-navigation-drawer>

    <v-main primary class="primary lighten-5 pl-14">
      <component v-bind:drawer.sync="drawer" :is="current"> </component>
    </v-main>
  </v-app>
</template>

<script>
import clustersLigandView from "./components/clustersLigandView";
import navDrawer from "./components/navDrawer";
import { mapState } from "vuex";
import infoOverlayComponent from "./components/infoOverlayComponent";

export default {
  name: "App",

  components: {
    clustersLigandView,
    navDrawer,
    infoOverlayComponent,
  },

  props: {
    header: String,
    text: String,
    iconColor: String,
  },

  data: () => ({
    // DOKU Overlay
    loading: false,
    //
    current: "clustersLigandView",
    navbarKey: true,
    drawer: true,
  }),
  mounted() {
    this.$store.dispatch("handleClusterData").then(() => {
      this.getFuncGroupData();
    });
    this.$store.dispatch("pollFGSMapData");
    this.$store.dispatch("pollGUIData");
    this.$store.dispatch("pollSelectionData");
    this.$store.dispatch("create_GUIcolorScaleSVGs");
  },
  computed: {
    ...mapState({
      showNavDrawer: (state) => state.showNavDrawer,
      infoBoxSegement: (state) => state.showInfoOverlay,
      showBlockUserInput: (state) => state.showBlockUserInput,
      infoOverlayContentHeader: (state) => state.infoOverlayContent[0],
      infoOverlayContentText: (state) => state.infoOverlayContent[1],
      isFullscreen: (state) => state.isFullscreen,
    }),
    getInfoOverlayState: {
      get() {
        return this.infoBoxSegement;
      },
    },
    getInfoOverlayContentText: {
      get() {
        return this.infoOverlayContent[1];
      },
    },
    getInfoOverlayContentHeader: {
      get() {
        return this.infoOverlayContent[0];
      },
    },
    isPinned: {
      get() {
        return !this.$store.state.showNavDrawer;
      },
      set(val) {
        console.log("val", val);
        this.$store.commit("setShowNavDrawer", !val);
      },
    },
    goFullscreen: {
      get() {
        var elem = document.getElementById("app");
        if (!this.isFullscreen) {
          if (document.exitFullscreen) {
            document.exitFullscreen();
          } else if (document.webkitExitFullscreen) {
            /* Safari */
            document.webkitExitFullscreen();
          } else if (elem.msExitFullscreen) {
            /* IE11 */
            document.msExitFullscreen();
          }
        } else {
          if (elem.requestFullscreen) {
            elem.requestFullscreen();
          } else if (elem.webkitRequestFullscreen) {
            /* Safari */
            elem.webkitRequestFullscreen();
          } else if (elem.msRequestFullscreen) {
            /* IE11 */
            elem.msRequestFullscreen();
          }
        }
        return this.isFullscreen;
      },
      set(val) {
        console.log("fullscreen", val);
        this.$store.commit("setIsFullscreen", val);
        var elem = document.getElementById("app");
        if (!val) {
          if (document.exitFullscreen) {
            document.exitFullscreen();
          } else if (document.webkitExitFullscreen) {
            /* Safari */
            document.webkitExitFullscreen();
          } else if (elem.msExitFullscreen) {
            /* IE11 */
            document.msExitFullscreen();
          }
        } else {
          if (elem.requestFullscreen) {
            elem.requestFullscreen();
          } else if (elem.webkitRequestFullscreen) {
            /* Safari */
            elem.webkitRequestFullscreen();
          } else if (elem.msRequestFullscreen) {
            /* IE11 */
            elem.msRequestFullscreen();
          }
        }
      },
    },
  },
  // methods of this vue component access via this.<funcName>
  methods: {
    clickAction(val) {
      this.$store.commit("setInfoOverlayState", val);
    },
    getFuncGroupData() {
      //console.log("parsing getFuncGroupData");
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

      fetch("/getFuncGroupData", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
      })
        .then((response) => response.json())
        .then((dataContent) => {
          this.$store.commit(
            "setFuncGroupData",
            dataContent.funcGroupData_recv
          );
          //console.log(dataContent)
        });
    },
  },
};
</script>
