import * as d3 from "d3";
import { rgb } from "d3";
//import * as d3collection from "d3-collection"
//import clusters from '@/components/clusters';
import { store, clusterDataKeys as k } from "../../main.js";

// fonzt size of stacked bar/ box plot labels
var fontSize = "11.5pt";

let onClickNoStake = function () {
  let clusterID = "All_withoutNoise";
  if (store.state.showNoise == true) {
    clusterID = "All";
  }
  store.commit("setClusterID", clusterID);
  console.log("set clusterID to", clusterID);
};

let selectCluster = function (elem, event, d, callee) {
  event.stopPropagation();
  if (d == null) {
    onClickNoStake();
    let clusterID = "All";
    highlightCluster(clusterID);
  } else if (d.data != null && callee == "bar") {
    let clusterID = d.data["id"];

    store.commit("setClusterID", clusterID);
    // uncomment code block if you want to select not only the cluster (+ ligand; + model)

    let lgdData = store.state.presentClustData;
    let id = d3.select(elem.parentNode).datum().key;
    store.commit("setLigandModelID", [
      id,
      [lgdData[id]["clusters"][d.data["id"]]["models"][0]],
    ]);
    highlightCluster(clusterID);
  } else {
    let clusterID = parseInt(d[0]);
    store.commit("setClusterID", clusterID);
    highlightCluster(clusterID);
  }
};

export function prepareIncomingData(dataset) {
  // Check if data available
  let preparedData = {};
  if (Object.keys(dataset).length === 0 && dataset.constructor === Object) {
    console.log("No data available. Please send data first.");
  } else {
    preparedData.mdlCount = dataset[k.modelCounts];
    preparedData.lgdIDs = dataset[k.ligandIDs];
    preparedData.lgdNames = dataset[k.zincNames];
    preparedData.assignedClust = dataset[k.clusterAssignment];
    preparedData.energies = dataset[k.energies];
    preparedData.mwt = dataset[k.mwt];
    preparedData.lipo = dataset[k.logp];
    preparedData.sp3 = dataset[k.fractioncsp3];
    preparedData.clusterSizes = dataset[k.clusterSizes];
    preparedData.mwtIdxSorted = Object.entries({ ...preparedData.mwt })
      .sort(([, a], [, b]) => a - b)
      .map((e) => e[0]);
    preparedData.lipoIdxSorted = Object.entries({ ...preparedData.lipo })
      .sort(([, a], [, b]) => a - b)
      .map((e) => e[0]);
    preparedData.sp3IdxSorted = Object.entries({ ...preparedData.sp3 })
      .sort(([, a], [, b]) => a - b)
      .map((e) => e[0]);

    preparedData.lgdNameID = {};
    preparedData.lgdEnergies = buildLgdInfoDict(dataset, k.energies);
    preparedData.hbonds = buildLgdInfoDict(dataset, k.hbonds);
    preparedData.halogenBonds = buildLgdInfoDict(dataset, k.halogenBonds);
    preparedData.hydrophobicInteractions = buildLgdInfoDict(
      dataset,
      k.hydrophobicInteractions
    );
    preparedData.metalComplexes = buildLgdInfoDict(dataset, k.metalComplexes);
    preparedData.piCationInteractions = buildLgdInfoDict(
      dataset,
      k.piCationInteractions
    );
    preparedData.piStacks = buildLgdInfoDict(dataset, k.piStacks);
    preparedData.saltBridges = buildLgdInfoDict(dataset, k.saltBridges);
    preparedData.efficiency = buildLgdInfoDict(dataset, k.efficiency);

    // build name/ID dict, needed for labelling the small barcharts
    for (let i = 0; i < preparedData.lgdNames.length; i++) {
      preparedData.lgdNameID[preparedData.lgdIDs[i]] = preparedData.lgdNames[i];
    }

    preparedData.uniqueClustIDs = [...new Set(preparedData.assignedClust)]; // gets all available clusterIDs

    let lgdData = {}; //contains names, clusters and models/nrgs contained in each cluster
    for (let i = 0; i < preparedData.lgdNames.length; i++) {
      lgdData[preparedData.lgdIDs[i]] = {};
      lgdData[preparedData.lgdIDs[i]]["name"] = preparedData.lgdNames[i];
      lgdData[preparedData.lgdIDs[i]]["mwt"] = preparedData.mwt[i];
      lgdData[preparedData.lgdIDs[i]]["lipo"] = preparedData.lipo[i];
      lgdData[preparedData.lgdIDs[i]]["sp3"] = preparedData.sp3[i];
      lgdData[preparedData.lgdIDs[i]]["clusters"] = {};
      lgdData[preparedData.lgdIDs[i]]["presentClustCnt"] = 0;
      lgdData[preparedData.lgdIDs[i]]["totalClustAvg"] = [];
      for (let c in preparedData.uniqueClustIDs) {
        lgdData[preparedData.lgdIDs[i]]["clusters"][
          preparedData.uniqueClustIDs[c]
        ] = { models: [], nrgs: [] };
      }
      lgdData[preparedData.lgdIDs[i]]["clusters"]["All"] = {
        models: [],
        nrgs: [],
      };
      lgdData[preparedData.lgdIDs[i]]["clusters"]["All_withoutNoise"] = {
        models: [],
        nrgs: [],
      };
    }

    // set lgdData in store as presentClustData for updating the clusters in ligand data table (see updateLgdTableData)
    store.commit("setPresentClustData", lgdData);
    store.commit("setClusterDataPrepared", preparedData);
  }
}
export function highlightCluster() {
  let color = "green";
  let heightHighlight = 0;
  let display = "block";
  let stackingDat = store.state.stackingDat;
  let highlightBar = store.state.highlightBar;
  let barScaleX = store.state.barScaleX;
  let barScaleY = store.state.barScaleY;
  let clusterID = store.state.clusterID;
  if (barScaleX == null) {
    return;
  }
  if (clusterID == "All" || clusterID == "All_withoutNoise") {
    display = "none";
  } else {
    let currentClusterData = stackingDat.filter((elem) => {
      return elem.id == clusterID;
    })[0];
    Object.keys(currentClusterData).forEach((key) => {
      if (key != "id" && key != "mean" && key != "median") {
        heightHighlight += currentClusterData[key];
        //console.log("Height Highlight", heightHighlight)
      }
    });
  }
  if (highlightBar) {
    highlightBar = d3
      .select("#highlightRect")
      .attr("x", barScaleX(clusterID) - 4)
      .attr("y", barScaleY(heightHighlight) - 4)
      .attr("width", barScaleX.bandwidth() + 8)
      .attr("height", barScaleY(0) - barScaleY(heightHighlight) + 8)
      .style("fill", "none")
      .style("stroke-width", "8")
      .style("stroke", color)
      .style("display", display);
  } else {
    highlightBar = d3
      .select("#highlightContainer")
      .append("rect")
      .attr("id", "highlightRect")
      .attr("x", barScaleX(clusterID) - 4)
      .attr("y", barScaleY(heightHighlight) - 4)
      .attr("width", barScaleX.bandwidth() + 8)
      .attr("height", barScaleY(0) - barScaleY(heightHighlight) + 8)
      .style("fill", "none")
      .style("stroke-width", "8")
      .style("stroke", color)
      .style("display", display);
  }
  store.commit("setHighlightBar", highlightBar);
}

export function buildLgdInfoDict(dataset, key) {
  let mdlCount = dataset[k.modelCounts];
  let info = dataset[key];
  let lgdIDs = dataset[k.ligandIDs];
  let lgdInfoDict = {};

  let indexCount = 0; // counter to get values

  for (let id in lgdIDs) {
    //assumes ligand id always starts with 0
    lgdInfoDict[id] = {};
    for (let mdl = 0; mdl < mdlCount[id]; mdl++) {
      // assumes models always start with 0
      lgdInfoDict[id][mdl] = {};
      lgdInfoDict[id][mdl] = info[indexCount];
      indexCount++;
    }
  }
  return lgdInfoDict;
}

function calculateQuantileScales(colors) {
  let mwt_q5 = d3.quantile(store.state.clusterData_prepared.mwt, 0.05);
  let mwt_q95 = d3.quantile(store.state.clusterData_prepared.mwt, 0.95);

  let lipo_q5 = d3.quantile(store.state.clusterData_prepared.lipo, 0.05);
  let lipo_q95 = d3.quantile(store.state.clusterData_prepared.lipo, 0.95);

  let sp3_q5 = d3.quantile(store.state.clusterData_prepared.sp3, 0.05);
  let sp3_q95 = d3.quantile(store.state.clusterData_prepared.sp3, 0.95);

  let quantileMWT = d3
    .scaleQuantize()
    .domain([mwt_q5, mwt_q95]) // pass the whole dataset to a scaleQuantile’s domain
    .range(colors.mwt);

  let quantileLipo = d3
    .scaleQuantize()
    .domain([lipo_q5, lipo_q95]) // pass the whole dataset to a scaleQuantile’s domain
    .range(colors.lipo);
  let quantileSp3 = d3
    .scaleQuantize()
    .domain([sp3_q5, sp3_q95]) // pass the whole dataset to a scaleQuantile’s domain
    .range(colors.sp3);
  return { mwt: quantileMWT, lipo: quantileLipo, sp3: quantileSp3 };
}

export function generateStackingData() {
  let preparedData = store.state.clusterData_prepared;
  let lgdData = store.state.presentClustData;
  // build data structure for proper stacking of the data
  let stackingDat = [];

  // initialize array[clusterCnt][ligandCnt]
  for (let i = 0; i < preparedData.uniqueClustIDs.length; i++) {
    // initialize clusters
    stackingDat[i] = { id: preparedData.uniqueClustIDs[i] };
    for (let j = 0; j < preparedData.lgdIDs.length; j++) {
      // initialise count for assigned clusters for each ligand
      stackingDat[i][preparedData.lgdIDs[j]] = 0;
    }
  }

  // initialise nrgDict to calculate mean energies of each cluster for sorting
  let nrgDict = {};
  for (let i = 0; i < stackingDat.length; i++) {
    nrgDict[stackingDat[i]["id"]] = {
      nrg: [],
      mean: 0,
      median: 0,
    };
  }

  let indexCount = 0; // counter to get energy values
  for (let id in preparedData.lgdIDs) {
    //assumes ligand id always starts with 0
    for (let mdl = 0; mdl < preparedData.mdlCount[id]; mdl++) {
      // assumes models always start with 0
      for (let cluster = 0; cluster < stackingDat.length; cluster++) {
        if (
          preparedData.assignedClust[indexCount] === stackingDat[cluster]["id"]
        ) {
          stackingDat[cluster][id] += 1;
          let datEnergy = preparedData.energies[indexCount];
          if (datEnergy != "" && datEnergy != 0) {
            nrgDict[stackingDat[cluster]["id"]]["nrg"].push(
              preparedData.energies[indexCount]
            ); // assumes energies and assigned clusters always have same length
            let curClusterEntry = stackingDat[cluster]["id"];
            lgdData[id]["clusters"][curClusterEntry]["models"].push(mdl);
            lgdData[id]["clusters"][curClusterEntry]["nrgs"].push(
              preparedData.energies[indexCount].toFixed(1)
            );
            if (parseInt(curClusterEntry) !== -1) {
              lgdData[id]["clusters"]["All_withoutNoise"]["models"].push(mdl);
              lgdData[id]["clusters"]["All_withoutNoise"]["nrgs"].push(
                preparedData.energies[indexCount].toFixed(1)
              );
            }
            lgdData[id]["clusters"]["All"]["models"].push(mdl);
            lgdData[id]["clusters"]["All"]["nrgs"].push(
              preparedData.energies[indexCount].toFixed(1)
            );
          }
        }
      }
      indexCount++;
    }
  }
  //console.log("lgdData", lgdData);

  // presentClustCount, totalClustAvg, does not take "none"/-1 cluster into account
  for (let l in lgdData) {
    for (let clust in lgdData[l]["clusters"]) {
      if (
        typeof lgdData[l]["clusters"] != "undefined" &&
        Object.keys(lgdData[l]["clusters"]).find(
          (key) => lgdData[l]["clusters"][key] === lgdData[l]["clusters"][clust]
        ) != -1
      ) {
        if (lgdData[l]["clusters"][clust]["models"].length > 0) {
          lgdData[l]["presentClustCnt"] += 1;
          for (
            let i = 0;
            i < lgdData[l]["clusters"][clust]["nrgs"].length;
            i++
          ) {
            lgdData[l]["totalClustAvg"].push(
              lgdData[l]["clusters"][clust]["nrgs"][i]
            );
          }
        }
      }
    }
  }

  // calculate mean, median, quartiles etc. for each cluster
  // round to 1 decimal to match input values

  for (let i in Object.keys(nrgDict)) {
    nrgDict[Object.keys(nrgDict)[i]]["mean"] = d3.mean(
      nrgDict[Object.keys(nrgDict)[i]]["nrg"]
    );
    nrgDict[Object.keys(nrgDict)[i]]["median"] = d3.median(
      nrgDict[Object.keys(nrgDict)[i]]["nrg"]
    );

    // for boxplot

    let q1 = d3.quantile(
      nrgDict[Object.keys(nrgDict)[i]]["nrg"].sort(d3.ascending),
      0.25
    );
    let q3 = d3.quantile(
      nrgDict[Object.keys(nrgDict)[i]]["nrg"].sort(d3.ascending),
      0.75
    );
    let interQuantileRange = q3 - q1;
    let min = d3.min(nrgDict[Object.keys(nrgDict)[i]]["nrg"]);
    let max = d3.max(nrgDict[Object.keys(nrgDict)[i]]["nrg"]);

    nrgDict[Object.keys(nrgDict)[i]]["q1"] = q1;
    nrgDict[Object.keys(nrgDict)[i]]["q3"] = q3;
    nrgDict[Object.keys(nrgDict)[i]]["iqr"] = interQuantileRange;
    nrgDict[Object.keys(nrgDict)[i]]["min"] = min;
    nrgDict[Object.keys(nrgDict)[i]]["max"] = max;
    nrgDict[Object.keys(nrgDict)[i]]["id"] = Object.keys(nrgDict)[i]; // needed for boxplot nesting
  }

  // store mean and median energies in stacking data
  for (let i = 0; i < stackingDat.length; i++) {
    stackingDat[i]["mean"] = nrgDict[stackingDat[i]["id"]]["mean"];
    stackingDat[i]["median"] = nrgDict[stackingDat[i]["id"]]["median"];
  }

  store.commit("setPresentClustData", lgdData);
  store.commit("setStackingDat", stackingDat);
  store.commit("setNrgDict", nrgDict);
}
export function drawWrapper(barDiv) {
  console.log("drawWrapper");
  drawBarchart(barDiv);
  drawBoxplot();
}

function drawBarchart(div) {
  console.log("draw the barchart");
  let coloring = store.state.stackedBar_coloring;
  //console.log("coloring",coloring);
  let sorting = store.state.sorting;
  let preparedData = store.state.clusterData_prepared;
  let lgdData = store.state.presentClustData;
  let stackingDat = store.state.stackingDat;
  let nrgDict = store.state.nrgDict;
  let showNoise = store.state.showNoise;
  let barScaleX = store.state.barScaleX;
  let barScaleY = store.state.barScaleY;
  store.commit("setHighlightBar", null);
  let maxVal = 0; // maxVal for y-axis

  for (let i in Object.keys(nrgDict)) {
    if (nrgDict[Object.keys(nrgDict)[i]]["nrg"].length > maxVal) {
      if (!showNoise && Object.keys(nrgDict)[i] === "-1");
      else {
        maxVal = nrgDict[Object.keys(nrgDict)[i]]["nrg"].length;
      }
    }
  }

  maxVal *= 1.06;

  // remove noise
  if (!showNoise) {
    stackingDat = stackingDat.filter((item) => item["id"] != -1);
  }

  // sort by mean energy
  let noiseCnt = preparedData.assignedClust.filter((x) => x == -1).length;
  if (sorting == "score") {
    stackingDat = stackingDat.sort(function (a, b) {
      console.log("sorting js ", a);
      return a.mean - b.mean;
    });
  } else if (sorting == "ligands") {
    stackingDat = stackingDat.sort(function (a, b) {
      let bCnt = b.id == "-1" ? noiseCnt : preparedData.clusterSizes[b.id];
      let aCnt = a.id == "-1" ? noiseCnt : preparedData.clusterSizes[a.id];
      return bCnt - aCnt;
    });
  } else if (sorting == "clusters") {
    stackingDat = stackingDat.sort(function (a, b) {
      return a.id - b.id;
    });
  }

  //console.log("coloring", coloring)
  let colors = {
    mwt: ["#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"],
    lipo: ["#ffffcc", "#c2e699", "#78c679", "#31a354", "#006837"],
    sp3: ["#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494"],
  };
  let quantileScales = calculateQuantileScales(colors);
  d3.select(div).selectAll("svg").remove();

  // Define width and height of the svg
  let margin = { top: 10, right: 30, bottom: 40, left: 60 };
  let box = document.querySelector("#wrapper").getBoundingClientRect();
  let svgWidth = box.width;
  let svgHeight = (box.height - (margin.top + margin.bottom)) * 0.65;

  // Append the svg object
  let svg_stackedBar = d3
    .select(div)
    .on("click", function (event, d) {
      selectCluster(this, event, d, "bar");
    })
    .classed("svgContainer", true)
    .append("svg")
    .attr("preserveAspectRatio", "xMinYMin meet")
    .attr(
      "viewBox",
      `0 0 ${svgWidth + margin.left + margin.right} ${
        svgHeight + margin.top + margin.bottom
      }`
    )
    .classed("resizeableSVG", true)
    //.attr("width", svgWidth + margin.left + margin.right)
    //.attr("height", svgHeight + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  d3.select(div)
    .append("div")
    .classed("chartTooltip", true)
    .attr("id", "tooltipClusterBar");
  // Append title
  //let title = "Assigned Clusters";
  svg_stackedBar
    .append("text")
    .attr("class", "plotTitle")
    .attr("x", svgWidth / 2 - margin.left / 2)
    .attr("y", 0 - margin.top / 2)
    .attr("text-anchor", "middle");
  //	.text(title);
  svg_stackedBar.append("g").attr("id", "highlightContainer");
  // Groups to be plotted
  stackingDat.columns = [];

  for (let i = 0; i < preparedData.lgdIDs.length; i++) {
    // TODO: FLEX selection
    if (coloring == "mwt")
      stackingDat.columns.push(preparedData.mwtIdxSorted[i].toString());
    else if (coloring == "lipo")
      stackingDat.columns.push(preparedData.lipoIdxSorted[i].toString());
    else if (coloring == "sp3")
      stackingDat.columns.push(preparedData.sp3IdxSorted[i].toString());
    else stackingDat.columns.push(preparedData.lgdIDs[i].toString());
  }
  let keys = stackingDat.columns.slice(0);

  /// x axis
  let xDomain = stackingDat.map(function (d) {
    return d.id;
  }); // Labeling clusters
  // sort domain in ascending order (sort by Cluster ID)

  barScaleX = d3
    .scaleBand()
    .domain(xDomain)
    .range([0, svgWidth - margin.left])
    .padding(0.5);
  store.commit("setBarScaleX", barScaleX);
  // customise ticks
  let tickVals = xDomain.map((x) => (x === -1 ? "noise" : "" + x.toString()));

  let xAxis = d3
    .axisBottom(barScaleX)
    .tickFormat((d, i) => tickVals[i])
    .tickSize(10);

  svg_stackedBar
    .append("g")
    .attr("transform", "translate(0," + svgHeight + ")")
    .call(xAxis)
    .selectAll("text")
    //.style("text-anchor", "end")
    //.attr("dx", "-1.6em")
    .attr("dy", "1.43em")
    .style("font-size", fontSize);
  //.attr("transform", "rotate(-90)"); //-65

  // x-label
  svg_stackedBar
    .append("text")
    .attr("class", "xlabel")
    .attr("text-anchor", "end")
    .attr("x", svgWidth)
    .attr("y", +svgHeight + 35)
    .style("font-size", fontSize)
    .text("ClusterID");

  // y axis
  barScaleY = d3.scaleLinear().domain([0, maxVal]).nice().range([svgHeight, 0]);
  store.commit("setBarScaleY", barScaleY);

  svg_stackedBar
    .append("g")
    .call(d3.axisLeft(barScaleY))
    .style("font-size", fontSize);

  // y-label
  svg_stackedBar
    .append("text")
    .attr("class", "ylabel")
    .attr("text-anchor", "end")
    .attr("transform", "rotate(-90)")
    .attr("y", -margin.left + 15)
    .attr("x", 0)
    .style("font-size", fontSize)
    .text("Ligand-Poses");

  // Color range
  let color = d3
    .scaleOrdinal()
    .domain(keys)
    .range([
      "#bf68b8",
      "#74b959",
      "#5a3789",
      "#c59a3a",
      "#6d80d8",
      "#79863a",
      "#b84873",
      "#45bc8d",
      "#ba533e",
      "#36dee6",
    ]);

  // Stack data
  let stackedData = d3.stack().keys(keys)(stackingDat);
  //console.log("keys", keys);
  //console.log("stackingDat", stackingDat);
  //console.log("stackedData", stackedData);

  // Tooltip
  let tooltip = d3
    .select("#tooltipClusterBar")
    .style("position", "absolute")
    .style("visibility", "hidden")
    .style("left", margin.left + 10 + "px");

  // Highlight bars of ligand subgroups on mouseover
  let mouseOver_stackedBar = function () {
    d3.select(this).style("cursor", "pointer");
    let ligand = d3.select(this.parentNode).datum().key;
    // set all bars to opacity = 0.1
    d3.selectAll(".ligand").style("opacity", 0.1);
    // set hovered bars/ligand back to opacity=1.0
    d3.selectAll(".ligand_" + ligand)
      .style("opacity", 1)
      .transition()
      .duration(0);
  };
  // All areas visible on mouseOut
  let mouseOut_stackedBar = function () {
    tooltip.style("visibility", "hidden");
    d3.selectAll(".ligand").style("opacity", 1);
  };

  // Changes tooltip text to be displayed
  let mouseMove = function (event, d) {
    let id = d3.select(this.parentNode).datum().key;
    let value = d.data[id];
    let modelIDs = [];

    lgdData[id]["clusters"][d.data["id"]]["models"].forEach((elem) => {
      modelIDs.push(parseInt(elem) + 1);
    });

    tooltip
      .style("visibility", "visible")
      .html(
        "<strong>Ligand " +
          id +
          ": </strong>#" +
          value +
          "<br>" +
          lgdData[id]["name"] +
          "<br><strong>Poses:</strong> " +
          modelIDs +
          "<br><strong>Energies:</strong> " +
          lgdData[id]["clusters"][d.data["id"]]["nrgs"] +
          "<br><strong>Molecular Weight:</strong> " +
          lgdData[id]["mwt"].toFixed(3) +
          "<br><strong>logP:</strong> " +
          lgdData[id]["lipo"].toFixed(3) +
          "<br><strong>Fraction sp<sup>3</sup>:</strong> " +
          lgdData[id]["sp3"].toFixed(3)
      )
      .style("left", event.layerX + 10 + "px")
      .style("top", event.layerY - 50 + "px");
  };

  // Plot areas
  svg_stackedBar
    .append("g")
    .attr("id", "mainGBar")
    .selectAll("g")
    .data(stackedData)
    .enter()
    .append("g")
    .style("shape-rendering", "crispEdges")
    .style("fill", function (d) {
      //TODO : flex color
      if (coloring == "mwt")
        return quantileScales.mwt([preparedData.mwt[d.key]]);
      else if (coloring == "lipo")
        return quantileScales.lipo([preparedData.lipo[d.key]]);
      else if (coloring == "sp3")
        return quantileScales.sp3([preparedData.sp3[d.key]]);
      else return color(d.key);
      //return color(d.key);
    })
    // add 2 classes to object
    .attr("class", function (d) {
      return "ligand ligand_" + d.key;
    })
    .selectAll("rect")
    .data(function (d) {
      return d;
    })
    .enter()
    .append("rect")
    .attr("x", function (d) {
      return barScaleX(d.data.id);
    })
    .attr("y", function (d) {
      return barScaleY(d[1]) - 1;
    })
    .attr("height", function (d) {
      return barScaleY(d[0]) - barScaleY(d[1]);
    })
    .attr("width", barScaleX.bandwidth)
    .on("mouseover", mouseOver_stackedBar)
    .on("mousemove", mouseMove)
    .on("click", function (event, d) {
      selectCluster(this, event, d, "bar");
    })
    .on("mouseleave", mouseOut_stackedBar);

  // plot highlight

  // legend adapted from: https://stackoverflow.com/questions/51520596/spread-d3-js-legend-on-two-columns
  //document.getElementById("legend").innerHTML = "";
  d3.select("#legend").selectAll("svg").remove();
  //console.log("drawSVGlegend");
  let svg2 = d3
    .select("#legend")
    .append("svg")
    .attr("width", 200)
    .attr("height", 20 * (preparedData.lgdIDs.length + 1))
    .append("g")
    .attr("transform", "translate(" + 0 + "," + 10 + ")");

  let legend = svg2
    .selectAll("legend")
    .data(stackedData)
    .enter()
    .append("g")
    .attr("transform", function (d, i) {
      return "translate(" + 0 + "," + Math.floor(i) * 20 + ")";
    });

  legend
    .append("text")
    .attr("x", 18)
    .attr("y", 10)
    .style("font-family", "Arial")
    .style("font-size", "10pt")
    .text(function (d) {
      return preparedData.lgdNameID[d.key];
    })
    .on("mouseover", mouseOver_stackedBar)
    .on("mouseleave", mouseOut_stackedBar);

  legend
    .append("rect")
    .attr("width", 10)
    .attr("height", 10)
    .attr("fill", function (d) {
      if (coloring == "mwt")
        return quantileScales.mwt([preparedData.mwt[d.key]]);
      else if (coloring == "lipo")
        return quantileScales.lipo([preparedData.lipo[d.key]]);
      else if (coloring == "sp3")
        return quantileScales.sp3([preparedData.sp3[d.key]]);
      else return color(d.key);
    })
    .on("mouseover", mouseOver_stackedBar)
    //.on("mousemove", showLgdName)
    .on("mouseleave", mouseOut_stackedBar);

  let colorModeLabel = function () {
    if (coloring === "lipo") {
      return "logP";
    } else if (coloring === "sp3") {
      return "Fsp3 ";
    } else if (coloring === "mwt") {
      return "Mwt ";
    } else {
      return "ID";
    }
  };
  let colorLegend = svg_stackedBar
    .append("text")
    .attr("x", 5)
    .attr("y", 12)
    .style("font-family", "Arial")
    .style("font-size", fontSize)
    .text("color mode: " + colorModeLabel());

  if (coloring != "id") {
    let spacing = 95;
    let offset_ColorModeLabel = 120;
    let quantilesScale = quantileScales[coloring].ticks();
    let valMinMax =
      coloring == "mwt"
        ? d3.extent(store.state.clusterData_prepared.mwt)
        : coloring == "lipo"
        ? d3.extent(store.state.clusterData_prepared.lipo)
        : d3.extent(store.state.clusterData_prepared.sp3);
    let toFixedNumber = 1;
    if (coloring == "mwt") {
      toFixedNumber = 0;
    }

    colorLegend = svg_stackedBar
      .selectAll()
      .data([
        {
          text: `< ${quantilesScale[0].toFixed(toFixedNumber)}`,
          color: colors[coloring][0],
        },
        {
          text: `${quantilesScale[0].toFixed(
            toFixedNumber
          )} – ${quantilesScale[1].toFixed(toFixedNumber)}`,
          color: colors[coloring][1],
        },
        {
          text: `${quantilesScale[1].toFixed(
            toFixedNumber
          )} – ${quantilesScale[2].toFixed(toFixedNumber)}`,
          color: colors[coloring][2],
        },
        {
          text: `${quantilesScale[2].toFixed(
            toFixedNumber
          )} – ${quantilesScale[3].toFixed(toFixedNumber)}`,
          color: colors[coloring][3],
        },
        {
          text: `> ${quantilesScale[3].toFixed(toFixedNumber)}`,
          color: colors[coloring][4],
        },
        {
          text: `min: ${valMinMax[0].toFixed(
            toFixedNumber
          )}, max: ${valMinMax[1].toFixed(toFixedNumber)} `,
          color: rgb(255, 255, 255),
        },
        //${valMinMax[0].toFixed(toFixedNumber)}
      ])
      .enter()
      .append("g")
      .attr("transform", "translate(" + 5 + "," + 0 + ")");

    colorLegend
      .append("text")
      .attr("x", function (d, i) {
        if (i == 0) {
          return 20 + offset_ColorModeLabel + 20;
        }
        return i * spacing + 20 + offset_ColorModeLabel;
      })
      .attr("y", 12)
      .style("font-family", "Arial")
      .style("font-size", fontSize)
      .text(function (d) {
        return d.text;
      });

    colorLegend
      .append("rect")
      .attr("x", function (d, i) {
        if (i == 0) {
          return 20 + offset_ColorModeLabel;
        }
        return i * spacing + offset_ColorModeLabel;
      })
      .attr("y", 0)
      .attr("width", 16)
      .attr("height", 16)
      .attr("fill", function (d) {
        return d.color;
      });
  }
  if (store.state.clusterID) {
    highlightCluster(store.state.clusterID);
  }
}

// adapted from https://www.d3-graph-gallery.com/graph/boxplot_several_groups.html
function drawBoxplot() {
  let nrgDict = store.state.nrgDict;
  let margin = { top: 10, right: 30, bottom: 40, left: 60 };
  let box = document.querySelector("#wrapper").getBoundingClientRect();
  let svgWidth = box.width;
  let svgHeight = (box.height - (margin.top + margin.bottom)) * 0.35;
  let showNoise = store.state.showNoise;
  let scaleX = store.state.barScaleX;
  let normalBoxCol = "lightblue";
  let boxHighlightCol = "royalblue";

  // Append the svg object
  let svg = d3
    .select("#assigned")
    .classed("svgContainer", true)
    .append("svg")
    .attr("preserveAspectRatio", "xMinYMin meet")
    .attr(
      "viewBox",
      `0 0 ${svgWidth + margin.left + margin.right} ${
        svgHeight + margin.top + margin.bottom
      }`
    )
    .classed("resizeableSVG2", true)
    //.attr("width", svgWidth + margin.left + margin.right)
    //.attr("height", svgHeight + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  d3.select("#assigned")
    .append("div")
    .classed("chartTooltip", true)
    .attr("id", "tooltipClusterBox");

  let dataReady = [];
  for (let i in Object.keys(nrgDict)) {
    dataReady.push(nrgDict[Object.keys(nrgDict)[i]]);
  }

  if (!showNoise) {
    dataReady = dataReady.filter((item) => item["id"] != "-1");
  }

  let dataGroup = d3.groups(dataReady, (d) => parseInt(d.id));
  let maxE = d3.max(
    dataReady.map(function (d) {
      return d["max"];
    })
  );
  // to get a bit space between box and axis

  let minE = d3.min(
    dataReady.map(function (d) {
      return d["min"];
    })
  );
  let rangePercentage = (minE - maxE) * 0.05;
  minE += rangePercentage * 2;
  maxE -= rangePercentage;

  //console.log("dataReady", dataReady);

  // generate empty ticks

  let tickVals = scaleX.domain().map(() => "");

  let xAxis = d3
    .axisTop(scaleX)
    .tickFormat((d, i) => tickVals[i])
    .tickSize(10);

  svg.append("g").call(xAxis);

  // y axis
  let scaleY = d3
    .scaleLinear()
    .domain([minE, maxE])
    .nice()
    .range([svgHeight, 0]);

  svg.append("g").call(d3.axisLeft(scaleY)).style("font-size", fontSize);

  // y-label
  svg
    .append("text")
    .attr("class", "ylabel")
    .attr("text-anchor", "start")
    .attr("transform", "rotate(-90)")
    .attr("y", -margin.left + 15)
    .attr("x", -svgHeight)
    .style("font-size", fontSize)
    .text("Score / Energy");

  // Tooltip
  let tooltip = d3
    .select("#tooltipClusterBox")
    .style("position", "absolute")
    .style("visibility", "hidden")
    .style("left", margin.left + 10 + "px");

  let mouseOver = function () {
    d3.select(this)
      .style("cursor", "pointer")
      .style("fill", boxHighlightCol)
      .style("opacity", 1)
      .style("stroke", "grey")
      .style("stroke-width", 2);
  };

  // Changes tooltip text to be displayed
  let mouseMove = function (event, d) {
    let mean = d[1][0].mean.toFixed(2);
    let median = d[1][0].median.toFixed(2);
    let q1 = d[1][0].q1.toFixed(1);
    let q3 = d[1][0].q3.toFixed(1);
    let maximum = d[1][0].min.toFixed(1); //most negative energy
    let minimum = d[1][0].max.toFixed(1); //least negative energy

    tooltip
      .style("visibility", "visible")
      .html(
        "<strong> Mean: </strong>" +
          mean +
          "<br>" +
          "<strong> Median: </strong>" +
          median +
          "<br>" +
          "<strong> Min: </strong>" +
          minimum +
          "<br>" +
          "<strong> Max: </strong>" +
          maximum +
          "<br>" +
          "<strong> Q1: </strong>" +
          q1 +
          "<br>" +
          "<strong> Q3: </strong>" +
          q3
      )
      .style("left", event.layerX + 10 + "px")
      .style("top", event.layerY - 100 + "px");
  };

  let mouseOut = function () {
    tooltip.style("visibility", "hidden");
    d3.select(this)
      .style("stroke", "black")
      .style("stroke-width", 1)
      .style("fill", normalBoxCol);
    //.style("opacity", 0.5);
  };

  let bw = scaleX.bandwidth();

  // vertical lines
  svg
    .selectAll("vertLines")
    .data(dataGroup)
    .enter()
    .append("line")
    .attr("x1", function (d) {
      return scaleX(d[0]) + bw / 2;
    })
    .attr("x2", function (d) {
      return scaleX(d[0]) + bw / 2;
    })
    .attr("y1", function (d) {
      return scaleY(d[1][0].min);
    })
    .attr("y2", function (d) {
      return scaleY(d[1][0].max);
    })
    .attr("stroke", "black")
    .style("width", 80);

  // boxes
  svg
    .selectAll("boxes")
    .data(dataGroup)
    .enter()
    .append("rect")
    .attr("x", function (d) {
      return scaleX(d[0]);
    })
    .attr("y", function (d) {
      return scaleY(d[1][0].q3);
    })
    .attr("height", function (d) {
      return scaleY(d[1][0].q1) - scaleY(d[1][0].q3);
    })
    .attr("width", bw)
    .attr("stroke", "black")
    .style("stroke-width", 1)
    .style("fill", normalBoxCol)
    //.style("opacity", 0.5)
    .on("mousemove", mouseMove)
    .on("mouseover", mouseOver)
    .on("mouseout", mouseOut)
    .on("click", function (event, d) {
      selectCluster(this, event, d, "box");
    });

  // median
  svg
    .selectAll("medianLines")
    .data(dataGroup)
    .enter()
    .append("line")
    .attr("x1", function (d) {
      return scaleX(d[0]);
    })
    .attr("x2", function (d) {
      return scaleX(d[0]) + bw;
    })
    .attr("y1", function (d) {
      return scaleY(d[1][0].median);
    })
    .attr("y2", function (d) {
      return scaleY(d[1][0].median);
    })
    .attr("stroke", "red")
    .style("width", 100);

  // mean
  svg
    .selectAll("meanLines")
    .data(dataGroup)
    .enter()
    .append("line")
    .attr("x1", function (d) {
      return scaleX(d[0]);
    })
    .attr("x2", function (d) {
      return scaleX(d[0]) + bw;
    })
    .attr("y1", function (d) {
      return scaleY(d[1][0].mean);
    })
    .attr("y2", function (d) {
      return scaleY(d[1][0].mean);
    })
    .attr("stroke", "blue")
    .style("width", 100);

  // whiskers max
  svg
    .selectAll("whiskersMax")
    .data(dataGroup)
    .enter()
    .append("line")
    .attr("x1", function (d) {
      return scaleX(d[0]) - bw / 3 + bw / 2;
    })
    .attr("x2", function (d) {
      return scaleX(d[0]) + bw / 3 + bw / 2;
    })
    .attr("y1", function (d) {
      return scaleY(d[1][0].max);
    })
    .attr("y2", function (d) {
      return scaleY(d[1][0].max);
    })
    .attr("stroke", "black")
    .style("width", 80);

  // whiskers min
  svg
    .selectAll("whiskersMin")
    .data(dataGroup)
    .enter()
    .append("line")
    .attr("x1", function (d) {
      return scaleX(d[0]) - bw / 3 + bw / 2;
    })
    .attr("x2", function (d) {
      return scaleX(d[0]) + bw / 3 + bw / 2;
    })
    .attr("y1", function (d) {
      return scaleY(d[1][0].min);
    })
    .attr("y2", function (d) {
      return scaleY(d[1][0].min);
    })
    .attr("stroke", "black")
    .style("width", 80);

  // legend
  let legend = svg
    .selectAll()
    .data(["Median", "Mean"])
    .enter()
    .append("g")
    .attr("transform", "translate(" + 0 + "," + (svgHeight - 10) + ")");
  //.attr("transform", "translate(" + (svgWidth-margin.left+10) + "," + (10) + ")");

  legend
    .append("text")
    .attr("y", 12)
    .attr("x", function (d) {
      return d === "Median" ? 18 : 86;
    })
    .style("font-family", "Arial")
    .style("font-size", fontSize)
    .text(function (d) {
      return d;
    });

  legend
    .append("rect")
    .attr("y", 2)
    .attr("x", function (d) {
      return d === "Median" ? 5 : 74;
    })
    .attr("width", 10)
    .attr("height", 10)
    .attr("fill", function (d) {
      return d === "Median" ? "red" : "blue";
    });
}

export function updateTable() {
  let lgdIDs = [];
  let clusterID = store.state.clusterID;
  let stackingData = store.state.stackingDat;
  let lgdData = store.state.presentClustData;
  let lgdEnergies = store.state.clusterData_prepared.lgdEnergies;

  for (let i = 0; i < stackingData.length; i++) {
    if (
      clusterID === stackingData[i].id ||
      clusterID === "All" ||
      clusterID === "All_withoutNoise"
    ) {
      for (let j = 0; j < Object.keys(stackingData[i]).length; j++) {
        let p = Array.from(Object.keys(stackingData[i]))[j].toString();
        if (
          stackingData[i][p] > 0 &&
          p != "id" &&
          p != "mean" &&
          p != "median"
        ) {
          lgdIDs.push(p);
        }
      }
    }
  }
  // delete duplicated entries
  lgdIDs = [...new Set(lgdIDs)];

  let tableData = [];
  let zincIDs = [];
  let mdlsInClust = {};

  let hbonds = store.state.clusterData_prepared.hbonds;
  let halogenBonds = store.state.clusterData_prepared.halogenBonds;
  let hydrophobicInteractions =
    store.state.clusterData_prepared.hydrophobicInteractions;
  let metalComplexes = store.state.clusterData_prepared.metalComplexes;
  let piCationInteractions =
    store.state.clusterData_prepared.piCationInteractions;
  let piStacks = store.state.clusterData_prepared.piStacks;
  let saltBridges = store.state.clusterData_prepared.saltBridges;

  // initialize  forceIsPresent
  let forceIsPresentPerLig = [];
  let forceIsPresent = [];
  for (const key of store.state.surfaceColoring_dropdownModesStore) {
    let item = key.text;
    forceIsPresent[item] = 0;
  }
  //console.log("forceIsPresent",forceIsPresent);

  // check if one of the models of a certain ligand contains a certain force;
  let cnt = 0;
  //console.log("lgdIDs",lgdIDs)
  for (let i of lgdIDs) {
    let curForPres = [...forceIsPresent];
    for (let j of lgdData[i]["clusters"][clusterID]["models"]) {
      //console.log("lig: ", i, "mdl: ", j);
      curForPres["HBonds"] = hbonds[i][j] > 0 ? "HBonds" : curForPres["HBonds"];
      curForPres["halogenBonds"] =
        halogenBonds[i][j] > 0 ? "halogenBonds" : curForPres["halogenBonds"];
      curForPres["hydrophobicInteractions"] =
        hydrophobicInteractions[i][j] > 0
          ? "hydrophobicInteractions"
          : curForPres["hydrophobicInteractions"];
      curForPres["metalComplexes"] =
        metalComplexes[i][j] > 0
          ? "metalComplexes"
          : curForPres["metalComplexes"];
      curForPres["piCationInteractions"] =
        piCationInteractions[i][j] > 0
          ? "piCationInteractions"
          : curForPres["piCationInteractions"];
      curForPres["piStacks"] =
        piStacks[i][j] > 0 ? "piStacks" : curForPres["piStacks"];
      curForPres["saltBridges"] =
        saltBridges[i][j] > 0 ? "saltBridges" : curForPres["saltBridges"];
      //console.log(curForPres);
    }

    forceIsPresentPerLig[cnt] = curForPres;
    cnt++;
  }
  //console.log(forceIsPresentPerLig);

  for (let i = 0; i < lgdIDs.length; i++) {
    let mdlAvgNrgs = [];
    mdlsInClust[lgdIDs[i]] =
      lgdData[lgdIDs[i]]["clusters"][clusterID]["models"];
    //console.log("mdlsInClust[lgdIDs[i]]", mdlsInClust[lgdIDs[i]]);
    for (let mdl in lgdEnergies[lgdIDs[i]]) {
      mdlAvgNrgs.push(lgdEnergies[lgdIDs[i]][mdl]);
    }
    let zincName = lgdData[lgdIDs[i]]["name"];
    zincIDs.push(zincName);
    let curClustAvg = d3.mean(
      lgdData[lgdIDs[i]]["clusters"][clusterID]["nrgs"]
    );
    if (typeof curClustAvg === "undefined") {
      curClustAvg = 0;
    }

    if (lgdData[lgdIDs[i]]["clusters"][clusterID]["models"].length != 0)
      tableData.push({
        lgd_id: lgdIDs[i],
        lgd_name: zincName,
        svg: zincName,
        mdl_cnt: lgdData[lgdIDs[i]]["clusters"][clusterID]["models"].length,
        local_min: Math.min(
          ...lgdData[lgdIDs[i]]["clusters"][clusterID]["nrgs"]
        ).toFixed(1),
        local_max: Math.max(
          ...lgdData[lgdIDs[i]]["clusters"][clusterID]["nrgs"]
        ).toFixed(1),
        avg_mdl: d3.mean(mdlAvgNrgs).toFixed(1), //average of all model energies
        clust_cnt: lgdData[lgdIDs[i]]["presentClustCnt"], // without "noise"
        curr_clust_avg: curClustAvg.toFixed(1), // mdlNrgs/mdlCnt in this cluster
        total_clust_avg:
          typeof d3.mean(lgdData[lgdIDs[i]]["totalClustAvg"]) === "undefined"
            ? "n/a"
            : d3.mean(lgdData[lgdIDs[i]]["totalClustAvg"]).toFixed(1), // mean of all models that are present in clusters (without "noise")
        forceIsPresentLig: forceIsPresentPerLig[i],
      });
  }

  store.commit("addSvgs", zincIDs);
  store.commit("setTableData", tableData);
  store.commit("setMdlInSelClust", mdlsInClust);
}

export function updateLgdTableData() {
  if (store.state.lgdTableTitle != "Ligand") {
    //console.log("updateLgdTableData")
    let tableHeaders2 = store.state.tableHeaders2;

    // ligandID, preparedData.lgdEnergies, lgdData, preparedData.hbonds, preparedData.efficiency
    let lgdId = store.state.ligandID;
    let lgdEnergies = store.state.lgdEnergies;
    let lgdData = store.state.presentClustData;
    let hbonds = store.state.clusterData_prepared.hbonds;
    let halogenBonds = store.state.clusterData_prepared.halogenBonds;
    let hydrophobicInteractions =
      store.state.clusterData_prepared.hydrophobicInteractions;
    let metalComplexes = store.state.clusterData_prepared.metalComplexes;
    let piCationInteractions =
      store.state.clusterData_prepared.piCationInteractions;
    let piStacks = store.state.clusterData_prepared.piStacks;
    let saltBridges = store.state.clusterData_prepared.saltBridges;
    let efficiency = store.state.clusterData_prepared.efficiency;

    let tableData2 = [];
    let presentClusters = [];
    let presentModels = [];
    let modelID = [];
    let modelIDs = [];

    for (modelID in lgdEnergies[lgdId]) {
      modelIDs.push(modelID);
      // iterate over clusters where a model of the ligandID is present
      for (let i = 0; i < Object.keys(lgdData[lgdId]["clusters"]).length; i++) {
        let clusterID = Object.keys(lgdData[lgdId]["clusters"])[i];
        // check which mdl of the ligand is present in the specific cluster
        //console.log(" lgdData[lgdId]['clusters'][clusterID]['models']", lgdData[lgdId]['clusters'][clusterID]['models'],"energy",modelID);
        if (
          lgdData[lgdId]["clusters"][clusterID]["models"].includes(
            parseInt(modelID)
          )
        ) {
          if (clusterID != "All" && clusterID != "All_withoutNoise") {
            presentClusters.push(clusterID);
          }
        }
      }

      //*** get models which contain selected func. group (if one is selected) ***//
      if (store.state.checkedFGSgroups != null) {
        if (store.state.checkB_fgsClust) {
          let s = store.state.curSelected_FGS_clusterData;

          for (let i = 0; i < s["fgsTypeIDs"].length; i++) {
            let short = s[s["fgsTypeIDs"][i]];
            //console.log("short", short);
            for (let i = 0; i < short["ligID"].length; i++) {
              let ligandID = short["ligID"][i];
              if (ligandID == lgdId) {
                //console.log("ligandID: ", ligandID);
                let mdlID = short["mdlID"][i];
                if (parseInt(modelID) == parseInt(mdlID)) {
                  if (presentModels[modelID] == null || -1) {
                    if (presentModels[modelID] == null)
                      presentModels.push(mdlID + 1);
                    if (presentModels[modelID] == -1)
                      presentModels[modelID] = mdlID + 1;
                  }
                } else {
                  if (presentModels[modelID] == null) presentModels.push(-1);
                }
              }
            }
          }
          //console.log("presentModels-fgsMapbased", presentModels);
        } else {
          store.state.checkedFGSgroups.forEach((element) => {
            if (
              store.state.clusterID == "All" ||
              store.state.clusterID == "All_withoutNoise"
            ) {
              for (const [clstrID] of Object.entries(
                store.state.clusterFuncData
              )) {
                //if(store.state.clusterID == "All_withoutNoise" && clstrID === -1) continue;
                let fgsCount = 0;
                let iter = 0;
                let value =
                  store.state.clusterFuncData[clstrID][element.fgsTypeID];
                if (value["funcGroupCnt"] > 0) {
                  fgsCount = value["funcGroupCnt"];
                  iter = value["lig_and_mdl_IDs"];
                  for (let i = 0; i < fgsCount; i++) {
                    let ligandID = iter[i * 2 + 0];
                    if (ligandID == lgdId) {
                      let mdlID = iter[i * 2 + 1];
                      if (parseInt(modelID) == parseInt(mdlID)) {
                        if (presentModels[modelID] == null || -1) {
                          if (presentModels[modelID] == null)
                            presentModels.push(mdlID + 1);
                          if (presentModels[modelID] == -1)
                            presentModels[modelID] = mdlID + 1;
                        }
                      } else {
                        if (presentModels[modelID] == null)
                          presentModels.push(-1);
                      }
                    }
                  }
                }
              }
            } else {
              //console.log("using clusterFuncData for filtering of model table")
              //filterZINC.push(...this.funcGroupData_inverted[element.fgsTypeID]);
              let fgsCount =
                store.state.clusterFuncData[store.state.clusterID][
                  element.fgsTypeID
                ]["funcGroupCnt"];
              let iter =
                store.state.clusterFuncData[store.state.clusterID][
                  element.fgsTypeID
                ]["lig_and_mdl_IDs"];
              for (let i = 0; i < fgsCount; i++) {
                let ligandID = iter[i * 2 + 0];
                if (ligandID == lgdId) {
                  let mdlID = iter[i * 2 + 1];
                  if (parseInt(modelID) == parseInt(mdlID)) {
                    if (presentModels[modelID] == null || -1) {
                      if (presentModels[modelID] == null)
                        presentModels.push(mdlID + 1);
                      if (presentModels[modelID] == -1)
                        presentModels[modelID] = mdlID + 1;
                    }
                  } else {
                    if (presentModels[modelID] == null) presentModels.push(-1);
                  }
                }
              }
            }
          });
        }
      }
    }
    //console.log("presentModels", presentModels);
    //console.log("presentClusters", presentClusters);

    for (let i = 0; i < modelIDs.length; i++) {
      let fingerPrint = [];
      //let no = "⚪";
      //let yes = "⚫";
      let no = false;
      let yes = true;
      fingerPrint[0] = hbonds[lgdId][modelIDs[i]] === 1 ? yes : no;
      fingerPrint[1] = halogenBonds[lgdId][modelIDs[i]] === 1 ? yes : no;
      fingerPrint[2] =
        hydrophobicInteractions[lgdId][modelIDs[i]] === 1 ? yes : no;
      fingerPrint[3] = metalComplexes[lgdId][modelIDs[i]] === 1 ? yes : no;
      fingerPrint[4] =
        piCationInteractions[lgdId][modelIDs[i]] === 1 ? yes : no;
      fingerPrint[5] = piStacks[lgdId][modelIDs[i]] === 1 ? yes : no;
      fingerPrint[6] = saltBridges[lgdId][modelIDs[i]] === 1 ? yes : no;
      tableData2.push({
        mdl_id: parseInt(parseInt(modelIDs[i]) + 1),
        mdl_nrg: lgdEnergies[lgdId][modelIDs[i]].toFixed(1),
        lgd_eff: efficiency[lgdId][modelIDs[i]].toFixed(1),
        presClust_id: presentClusters[i] == -1 ? "noise" : presentClusters[i],
        presModel_id: presentModels[i] == null ? -1 : presentModels[i],
        hbond: fingerPrint,
        /*
				halogenBonds: (halogenBonds[lgdId][modelIDs[i]] === 1 ? true : false),
				hydrophobicInteractions: (hydrophobicInteractions[lgdId][modelIDs[i]] === 1 ? true : false),
				metalComplexes: (metalComplexes[lgdId][modelIDs[i]] === 1 ? true : false),
				piCationInteractions: (piCationInteractions[lgdId][modelIDs[i]] === 1 ? true : false),
				piStacks: (piStacks[lgdId][modelIDs[i]] === 1 ? true : false),
				saltBridges: (saltBridges[lgdId][modelIDs[i]] === 1 ? true : false),
				*/

        //hbond: (hbonds[lgdId][modelIDs[i]] === 1 ? "✅" : "❌"),
      });
    }
    store.commit("setTableData2", tableData2);
    store.commit("setTableHeaders2", tableHeaders2);
  }
}

export function parseFGSMapData(dataset) {
  let curSelected_FGS_clusterData = [];
  curSelected_FGS_clusterData["fgsTypeIDs"] = [];
  let i = 0;
  //console.log("parsing FGSMapData");
  while (i < dataset.length) {
    let fgsTypeID = dataset[i];
    i++;
    let oneFGSTypeArray = [];
    curSelected_FGS_clusterData["fgsTypeIDs"].push(fgsTypeID);
    oneFGSTypeArray["ligID"] = [];
    oneFGSTypeArray["mdlID"] = [];
    while (dataset[i] != -1 && i < dataset.length) {
      oneFGSTypeArray["ligID"].push(dataset[i + 0]);
      oneFGSTypeArray["mdlID"].push(dataset[i + 1]);
      i += 2;
    }
    curSelected_FGS_clusterData[fgsTypeID] = oneFGSTypeArray;
    i++;
  }
  //console.log("curSelected_FGS_clusterData",curSelected_FGS_clusterData);
  return curSelected_FGS_clusterData;
}
