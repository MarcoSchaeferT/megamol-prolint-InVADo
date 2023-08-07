<template>
  <v-list nav dense>
    <v-expansion-panels>
      <v-expansion-panel>
        <v-expansion-panel-header>
          <div class="font-weight-bold">Protein Data</div>
        </v-expansion-panel-header>
        <v-expansion-panel-content>
          <v-checkbox
            v-model="GUI"
            label="cartoon represeantation"
            value="isCartoonRenderer"
          ></v-checkbox>
          <v-checkbox
            v-model="GUI"
            label="binding site area (only)"
            value="pocketArea"
          ></v-checkbox>
          <v-checkbox
            v-model="GUI"
            label="transparency"
            value="transparency"
          ></v-checkbox>
          <v-row
            ><v-col cols="1"></v-col
            ><v-col>
              <v-checkbox
                v-model="GUI"
                label="opaque binding site"
                value="opaquePocket"
              ></v-checkbox> </v-col
          ></v-row>
          <v-checkbox
            v-model="GUI"
            label="protein residues (on/off)"
            value="residues"
          ></v-checkbox>
          <v-spacer></v-spacer>
          <v-row pa-0>
            <v-col cols="1"></v-col>
            <v-col>
              <v-checkbox
                v-model="GUI"
                label="protein backbone (on/off)"
                value="residuesBackbone"
              ></v-checkbox>
              <v-checkbox
                v-model="GUI"
                label="protein residue labels (on/off)"
                value="residuesLabels"
              ></v-checkbox>
              <v-checkbox
                v-model="GUI"
                label="protein residues (interaction only)"
                value="interactionResidues"
              ></v-checkbox>
            </v-col>
          </v-row>
          <v-checkbox
            v-model="GUI"
            label="chemical property"
            value="chemicalProperty"
          ></v-checkbox>
          <v-row
            ><v-col cols="1"> </v-col
            ><v-col>
              <v-select
                :items="proteinColoringModes"
                v-model="proteinColoringMode0"
                label="50% --- protein coloring mode 0"
                @click="setState"
                @input="setStateBck"
                outlined
              ></v-select>
              <v-select
                :items="proteinColoringModes"
                v-model="proteinColoringMode1"
                label="50% --- protein coloring mode 1"
                @click="setState"
                @input="setStateBck"
                outlined
              ></v-select> </v-col
          ></v-row>
        </v-expansion-panel-content>
      </v-expansion-panel>
      <v-expansion-panel>
        <v-expansion-panel-header>
          <div class="font-weight-bold">Ligand Data</div>
        </v-expansion-panel-header>
        <v-expansion-panel-content>
          <v-checkbox
            v-model="GUI"
            label="show interaction type bar chart (on/off)"
            value="barchartGylphe"
          ></v-checkbox>
          <v-checkbox
            v-model="GUI"
            label="show cluster dodecahedrons"
            value="clusterSpheres"
          ></v-checkbox>
          <div v-html="getGUIcolorScale('pocketCluster')"></div>
          <v-checkbox
            v-model="GUI"
            label="show functional group cluster"
            value="funcGroupCentroids"
          ></v-checkbox>
          <div v-html="getGUIcolorScale('funcGroupCluster')"></div>
          <v-checkbox
            v-model="GUI"
            label="show all cluster atoms"
            value="clusterAtoms"
          ></v-checkbox>
          <v-checkbox
            v-model="GUI"
            label="show noise binding poses"
            value="clusterNoise"
          ></v-checkbox>
          <!--<v-checkbox
            v-model="GUI"
            label="improve text readability"
            value="improveTextReadability"
          ></v-checkbox>
        -->
        </v-expansion-panel-content>
      </v-expansion-panel>
      <v-expansion-panel>
        <v-expansion-panel-header>
          <div class="font-weight-bold">Interaction Types:</div>
          Surface Coloring
        </v-expansion-panel-header>
        <v-expansion-panel-content>
          <v-checkbox
            v-model="GUI"
            label="interaction types (on/off)"
            value="interactionForce"
          ></v-checkbox>
          <v-select
            :items="surfaceColoring_dropdownModes"
            v-model="surfaceColoring_dropdown"
            label="surface coloring force"
            @click="setState"
            @input="setStateBck"
            outlined
            hide-details=""
            class="pa-0"
          >
            <template #selection="{ item }">
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
          <v-checkbox
            v-model="GUI"
            label="color by count of interaction type"
            value="interactionForceAmountColoring"
          ></v-checkbox>
          <div v-html="getGUIcolorScale('interactionForce')"></div>
        </v-expansion-panel-content>
      </v-expansion-panel>
      <v-expansion-panel>
        <v-expansion-panel-header>
          <div class="font-weight-bold">Interaction Types:</div>
          Sticks
        </v-expansion-panel-header>
        <v-expansion-panel-content>
          <v-checkbox v-model="GUI" label="hydrogenBonds" value="hydrogenBonds">
            <template v-slot:label>
              <span style="color: rgb(103.785, 49.98, 167.79)">&#11044;</span
              >&nbsp;hydrogen bonds
            </template>
          </v-checkbox>
          <v-checkbox v-model="GUI" label="hydrogenCones" value="hydrogenCones"
            ><template v-slot:label>
              <!--<span style="color: rgb(228, 198, 242)">&#11044;</span>-->
              <span style="color: rgb(142, 53, 155)">&#11044;</span
              >&nbsp;hydrogen cones
            </template>
          </v-checkbox>
          <v-checkbox v-model="GUI" label="halogenBonds" value="halogenBonds"
            ><template v-slot:label>
              <span style="color: rgb(55.08, 125.97, 184.11)">&#11044;</span
              >&nbsp;halogen bonds
            </template>
          </v-checkbox>
          <v-checkbox
            v-model="GUI"
            label="hydrophobicInteractions"
            value="hydrophobicInteractions"
            ><template v-slot:label>
              <span style="color: rgb(77.01, 174.93, 73.95)">&#11044;</span
              >&nbsp;hydrophobic interactions
            </template>
          </v-checkbox>
          <v-checkbox
            v-model="GUI"
            label="metalComplexes"
            value="metalComplexes"
            ><template v-slot:label>
              <span style="color: rgb(166, 86, 40)">&#11044;</span>&nbsp;metal
              complexes
            </template>
          </v-checkbox>
          <v-checkbox
            v-model="GUI"
            label="piCationInteractions"
            value="piCationInteractions"
            ><template v-slot:label>
              <span style="color: rgb(255, 126.99, 0)">&#11044;</span
              >&nbsp;&pi;-cation interactions
            </template>
          </v-checkbox>
          <v-checkbox v-model="GUI" label="piStacks" value="piStacks"
            ><template v-slot:label>
              <span style="color: rgb(245, 221, 25)">&#11044;</span
              >&nbsp;&pi;-Stacks
            </template>
          </v-checkbox>
          <v-checkbox v-model="GUI" label="saltBridges" value="saltBridges"
            ><template v-slot:label>
              <span style="color: rgb(226.95, 26.01, 28.05)">&#11044;</span
              >&nbsp;salt bridges
            </template>
          </v-checkbox>
        </v-expansion-panel-content>
      </v-expansion-panel>

      <v-expansion-panel>
        <v-expansion-panel-header>
          <div class="font-weight-bold">Clustering</div>
          Ligand Binding Poses
        </v-expansion-panel-header>
        <v-expansion-panel-content>
          <v-col class="d-flex align-center" cols="9">
            <v-text-field
              hide-details
              :value="modelClustering_eps"
              @change="
                (event) => {
                  modelClustering_eps = event;
                  clusteringON = true;
                }
              "
              label="search radius (Å)"
              class="mt-0 pt-0"
              type="number"
              style="width: 60px"
            ></v-text-field
          ></v-col>

          <v-col class="d-flex align-center" cols="9">
            <v-text-field
              hide-details
              :value="modelClustering_minBindEnergy"
              @change="
                (event) => {
                  modelClustering_minBindEnergy = event;
                  clusteringON = true;
                }
              "
              label="min binding score/energy"
              class="mt-0 pt-0"
              type="number"
              style="width: 60px"
            ></v-text-field
          ></v-col>

          <v-col class="d-flex align-center" cols="9">
            <v-text-field
              hide-details
              :value="modelClustering_minPts"
              @change="
                (event) => {
                  modelClustering_minPts = event;
                  clusteringON = true;
                }
              "
              label="min items in cluster"
              class="mt-0 pt-0"
              type="number"
              style="width: 60px"
            ></v-text-field
          ></v-col>

          <v-progress-circular
            v-if="clusteringON"
            indeterminate
            color="primary"
          ></v-progress-circular>
          <div style="color: darkgrey; text-size: 8pt">
            *change and press Enter
          </div>
        </v-expansion-panel-content>
      </v-expansion-panel>

      <v-expansion-panel>
        <v-expansion-panel-header>
          <div class="font-weight-bold">Clustering</div>
          Functional Groups
        </v-expansion-panel-header>
        <v-expansion-panel-content>
          <v-col class="d-flex align-center" cols="9">
            <v-text-field
              hide-details
              :value="searchRadius_fgc"
              @change="
                (event) => {
                  searchRadius_fgc = event;
                  clusteringON = true;
                }
              "
              label="search radius (Å)"
              class="mt-0 pt-0"
              type="number"
              style="width: 60px"
            ></v-text-field
          ></v-col>

          <v-col class="d-flex align-center" cols="9">
            <v-text-field
              hide-details
              v-model="minPTS_fgc"
              @change="
                () => {
                  clusteringON = true;
                }
              "
              label="min items in cluster"
              class="mt-0 pt-0"
              type="number"
              style="width: 60px"
            ></v-text-field
          ></v-col>
          <v-progress-circular
            v-if="clusteringON"
            indeterminate
            color="primary"
          ></v-progress-circular>
          <div style="color: darkgrey; text-size: 8pt">
            *change and press Enter
          </div>
        </v-expansion-panel-content>
      </v-expansion-panel>

      <v-expansion-panel>
        <v-expansion-panel-header>
          <div class="font-weight-bold">Radial Menu</div>
        </v-expansion-panel-header>
        <v-expansion-panel-content>
          <v-col class="d-flex align-center" cols="9">
            <v-text-field
              hide-details
              v-model="starMenuTextSize"
              label="menu text size"
              class="mt-0 pt-0"
              type="number"
              style="width: 60px"
            ></v-text-field
          ></v-col>

          <v-col class="d-flex align-center" cols="9">
            <v-text-field
              hide-details
              v-model="starMenuCircleCnt"
              label="radial menu circle count"
              class="mt-0 pt-0"
              type="number"
              style="width: 60px"
            ></v-text-field
          ></v-col>

          <v-col class="d-flex align-center" cols="9">
            <v-text-field
              hide-details
              v-model="starMenuSize"
              label="radial menu size"
              class="mt-0 pt-0"
              type="number"
              style="width: 60px"
            ></v-text-field
          ></v-col>

          <v-col class="d-flex align-center" cols="9">
            <v-text-field
              hide-details
              v-model="subSphereSizefactor"
              label="min circle size"
              class="mt-0 pt-0"
              type="number"
              style="width: 60px"
            ></v-text-field
          ></v-col>
          <v-select
            :items="starMenuSortingForces_dropdownModes"
            v-model="starMenuSorting_dropdown"
            label="sorting force"
            @click="setState"
            @input="setStateBck"
            outlined
          >
            <template #selection="{ item }">
              <span :style="'color:' + item.color">&#11044;</span>&nbsp;{{
                item.text
              }}
            </template>
            <template #item="{ item }">
              <span :style="'color:' + item.color">&#11044;</span>&nbsp;{{
                item.text
              }}
            </template></v-select
          >
        </v-expansion-panel-content>
      </v-expansion-panel>

      <v-expansion-panel>
        <v-expansion-panel-header>
          <div class="font-weight-bold">Clip Plane</div>
        </v-expansion-panel-header>
        <v-expansion-panel-content>
          <v-checkbox
            v-model="GUI"
            label="clip plane (on/off)"
            value="clipPlaneEnable"
          ></v-checkbox>
          <v-select
            hide-details
            :items="planeOrientations"
            v-model="planeOrientation"
            label="clip plane orientation"
            @click="setState"
            @input="setStateBck"
            outlined
          ></v-select>
          <v-subheader>clip plane percent</v-subheader>
          <v-col class="d-flex align-center" cols="12">
            0%
            <v-slider
              :value="clipPlaneDist"
              thumb-label="always"
              thumb-color="orange darken-3"
              thumb-size="20"
              hide-details
              :max="100"
              :min="0"
              @change="
                (event) => {
                  clipPlaneDist = event;
                }
              "
              class="mt-3 pt-0"
              type="number"
              style="width: 60px"
            ></v-slider
            >100%
          </v-col>
        </v-expansion-panel-content>
      </v-expansion-panel>
    </v-expansion-panels>
  </v-list>
</template>

<script>
import { mapState } from "vuex";

export default {
  // name of the component
  name: "navDrawer",
  //props: ['drawer'],
  components: {},

  // data section of the Vue component. Access via this.<varName> .
  data: () => ({
    // initial values set as in ModelClusterRenderer
    // checkboxes

    boxStore: null,
    clusteringON: false,
  }),

  computed: {
    ...mapState({
      //boxStore: state => state.GUIdata,
      /*
       * chemicalProperty:
       * clusterAtoms:
       * clusterNoise:
       * clusterSpheres:
       * funcGroupCentroids:
       * hydrogenBonds:
       * hydrogenCones:
       * improveTextReadability:
       * interactionForce:
       * interactionResidues:
       * minPTS_fgc:
       * modelClustering_eps:
       * modelClustering_minBindEnergy:
       * modelClustering_minPts:
       * opaquePocket:
       * pocketArea:
       * proteinColoringMode0:
       * proteinColoringMode1:
       * residues:
       * residuesLabels:
       * searchRadius_fgc:
       * transparency:
       */
      surfaceColoring_dropdownModesStore: (state) =>
        state.surfaceColoring_dropdownModesStore,
      starMenuSortingForces_dropdownModesStore: (state) =>
        state.starMenuSortingForces_dropdownModesStore,
      showNavDrawerVal: (state) => state.showNavDrawer,
      GUIcolorScaleSVGs: (state) => state.GUIcolorScaleSVGs,
    }),

    icon() {
      if (this.likesAllFruit) return "mdi-close-box";
      if (this.likesSomeFruit) return "mdi-minus-box";
      return "mdi-checkbox-blank-outline";
    },

    GUI: {
      get: function () {
        // eslint-disable-next-line vue/no-side-effects-in-computed-properties
        return this.$store.state.GUIdata["boolNames"];
      },
      set: function (val) {
        console.log("val", val);
        if (this.boxStore == null) {
          this.boxStore = Object.assign(this.$store.state.GUIdata);
        }
        if (this.boxStore["boolNames"].length > val.length) {
          let difference = this.boxStore["boolNames"].filter(
            (x) => !val.includes(x)
          );
          console.log("val= ", difference[0], ": false");
          this.postNavElem([difference[0], 0]);
        } else {
          if (val[val.length - 1] == "funcGroupCentroids") {
            this.$store.commit("setCheckB_fgsClust", true);
          }
          console.log("val= ", val[val.length - 1], ": true");
          this.postNavElem([val[val.length - 1], 1]);
        }
        console.log("val boxStore", this.boxStore["boolNames"]);
        let GUIdata = this.$store.state.GUIdata;
        GUIdata["boolNames"] = val;
        this.boxStore = GUIdata;
        this.$store.commit("setGUIdata", GUIdata);
      },
    },
    // floats
    modelClustering_eps: {
      get: function () {
        return this.roundFloat(
          this.checkGUIdata(
            this.$store.state.GUIdata["float"],
            "modelClustering_eps"
          )
        );
      },
      set: function (val) {
        this.postNavElem(["modelClustering_eps", val]);
        //let GUIdata = this.$store.state.GUIdata;
        //GUIdata["float"]["modelClustering_eps"] = val;
        //this.$store.commit("setGUIdata", GUIdata);
      },
    },

    modelClustering_minBindEnergy: {
      get: function () {
        return this.roundFloat(
          this.checkGUIdata(
            this.$store.state.GUIdata["float"],
            "modelClustering_minBindEnergy"
          )
        );
      },
      set: function (val) {
        this.postNavElem(["modelClustering_minBindEnergy", val]);
      },
    },
    searchRadius_fgc: {
      get: function () {
        return this.roundFloat(
          this.checkGUIdata(
            this.$store.state.GUIdata["float"],
            "searchRadius_fgc"
          )
        );
      },
      set: function (val) {
        this.postNavElem(["searchRadius_fgc", val]);
      },
    },
    clipPlaneDist: {
      get: function () {
        return this.checkGUIdata(
          this.$store.state.GUIdata["float"],
          "clipPlaneDist"
        );
      },
      set: function (val) {
        console.log("tzr", val);
        this.postNavElem(["clipPlaneDist", val]);
      },
    },
    //ints
    starMenuTextSize: {
      get: function () {
        return this.checkGUIdata(
          this.$store.state.GUIdata["int"],
          "starMenuTextSize"
        );
      },
      set: function (val) {
        this.postNavElem(["starMenuTextSize", val]);
      },
    },
    subSphereSizefactor: {
      get: function () {
        return this.checkGUIdata(
          this.$store.state.GUIdata["int"],
          "subSphereSizefactor"
        );
      },
      set: function (val) {
        this.postNavElem(["subSphereSizefactor", val]);
      },
    },
    modelClustering_minPts: {
      get: function () {
        return this.checkGUIdata(
          this.$store.state.GUIdata["int"],
          "modelClustering_minPts"
        );
      },
      set: function (val) {
        this.postNavElem(["modelClustering_minPts", val]);
      },
    },
    minPTS_fgc: {
      get: function () {
        return this.checkGUIdata(
          this.$store.state.GUIdata["int"],
          "minPTS_fgc"
        );
      },
      set: function (val) {
        this.postNavElem(["minPTS_fgc", val]);
      },
    },
    starMenuCircleCnt: {
      get: function () {
        return this.checkGUIdata(
          this.$store.state.GUIdata["int"],
          "starMenuCircleCnt"
        );
      },
      set: function (val) {
        this.postNavElem(["starMenuCircleCnt", val]);
      },
    },
    starMenuSize: {
      get: function () {
        return this.checkGUIdata(
          this.$store.state.GUIdata["int"],
          "starMenuSize"
        );
      },
      set: function (val) {
        this.postNavElem(["starMenuSize", val]);
      },
    },
    proteinColoringMode0: {
      get: function () {
        return this.checkGUIdata(
          this.$store.state.GUIdata["int"],
          "proteinColoringMode0"
        );
      },
      set: function (val) {
        this.postNavElem(["proteinColoringMode0", val]);
      },
    },
    proteinColoringMode1: {
      get: function () {
        return this.checkGUIdata(
          this.$store.state.GUIdata["int"],
          "proteinColoringMode1"
        );
      },
      set: function (val) {
        this.postNavElem(["proteinColoringMode1", val]);
      },
    },
    surfaceColoring_dropdown: {
      get: function () {
        return this.checkGUIdata(
          this.$store.state.GUIdata["int"],
          "surfaceColoring_dropdown"
        );
      },
      set: function (val) {
        this.postNavElem(["surfaceColoring_dropdown", val]);
      },
    },
    starMenuSorting_dropdown: {
      get: function () {
        return this.checkGUIdata(
          this.$store.state.GUIdata["int"],
          "starMenuSorting_dropdown"
        );
      },
      set: function (val) {
        this.postNavElem(["starMenuSorting_dropdown", val]);
      },
    },
    proteinColoringModes: {
      get: function () {
        return [
          { text: "ELEMENT", value: 0 },
          { text: "STRUCTURE", value: 1 },
          { text: "RAINBOW", value: 2 },
          { text: "BFACTOR", value: 3 },
          { text: "CHARGE", value: 4 },
          { text: "OCCUPANCY", value: 5 },
          { text: "CHAIN", value: 6 },
          { text: "MOLECULE", value: 7 },
          { text: "RESIDUE", value: 8 },
          // {text: "CHAINBOW"			, value: 9},
          { text: "AMINOACID", value: 10 },
          //{text: "VALUE"				, value: 11},
          //  {text: "CHAIN_ID"			, value: 12},
          // {text: "MOVEMENT"			, value: 13},
          { text: "BINDINGSITE", value: 14 },
          { text: "HYDROPHOBICITY", value: 15 },
          //{text: "PER_ATOM_FLOAT"	, value: 16},
          { text: "HEIGHTMAP_COL", value: 17 },
          // {text: "HEIGHTMAP_VAL"		, value: 18}
        ];
      },
    },
    surfaceColoring_dropdownModes: {
      get: function () {
        return this.surfaceColoring_dropdownModesStore;
      },
    },
    starMenuSortingForces_dropdownModes: {
      get: function () {
        return this.starMenuSortingForces_dropdownModesStore;
      },
    },
    planeOrientation: {
      get: function () {
        return this.checkGUIdata(
          this.$store.state.GUIdata["int"],
          "clipPlaneOrientation"
        );
      },
      set: function (val) {
        this.postNavElem(["clipPlaneOrientation", val]);
      },
    },
    planeOrientations: {
      get: function () {
        return [
          { text: "XY-plane", value: 0 },
          { text: "XZ-plane", value: 1 },
          { text: "YZ-plane", value: 2 },
        ];
      },
    },
  },

  // Setup a Listener for one or more variables of "data", i.e. slider1 references to data -> slider1
  watch: {
    // checkboxes
    chemicalProperty: function () {
      this.postNavElem(["chemicalProperty", this.chemicalProperty]);
    },
  },

  // functions to call on mount (after DOM etc. is built)
  mounted() {},

  // methods of this vue component access via this.<funcName>
  methods: {
    checkGUIdata(val, name) {
      if (val == undefined) {
        return 0;
      } else {
        return val[name];
      }
    },
    getGUIcolorScale(scaleName) {
      return this.GUIcolorScaleSVGs[scaleName];
    },
    roundFloat(val) {
      //val = Math.round(val * 10) / 10;
      val = Number(val.toFixed(1));
      return val;
    },
    postNavElem(val) {
      let GUI = this.$store.state.GUIdata;
      // convert bools to numbers
      if (typeof val[1] === "boolean") {
        val[1] === true ? (val[1] = 1) : (val[1] = 0);
        GUI["boolNames"][val[0]] = val[1];
      } else if (typeof val[1] === "number" && val[1] % 1 !== 0) {
        GUI["float"][val[0]] = val[1];
      } else if (typeof val[1] === "number" && val[1] % 1 === 0) {
        GUI["int"][val[0]] = val[1];
      }
      console.log("set: ", val[0], " val: ", val[1]);
      console.log("ALTT_:", GUI);
      this.$store.commit("setGUIdata", GUI);
      let is_showNavDrawer_alreadySetTrue = this.showNavDrawerVal;
      this.$store.commit("setShowNavDrawer", true);
      this.$store.commit("setBlockUserInput", true);
      fetch("/navigation", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({ exportData: val }),
      }).then(() => {
        if (!is_showNavDrawer_alreadySetTrue) {
          this.$store.commit("setShowNavDrawer", false);
        }
      });
    },

    setState() {
      console.log("fire state");
      this.$store.commit("setShowNavDrawer", true);
    },
    setStateBck() {
      console.log("fire state");
      this.$store.commit("setShowNavDrawer", false);
    },
  },
};
</script>
