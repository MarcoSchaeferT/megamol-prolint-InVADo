// partly from: https://www.d3-graph-gallery.com/graph/histogram_basic.html
import * as d3 from "d3";
import { store, clusterDataKeys as k } from "../../main.js";

/***
 * main histogram function
 * @param data
 * @param lgdNameID
 * @param lgdEnergies
 */

export function drawLigandHisto(data, binNum, ligandModel) {
  let margin = { top: 10, right: 20, bottom: 25, left: 40 };

  // Define width and height of the svg
  let box = document.querySelector("#ligand").getBoundingClientRect();
  let svgWidth = box.width - margin.left - margin.right;
  let svgHeight = 250 - margin.top - margin.bottom;

  // Check if data available
  if (Object.keys(data).length === 0 && data.constructor === Object) {
    console.log("No data available. Please send data first.");
  } else {
    console.log("mwt: ", data[k.mwt]);
    console.log("logp: ", data[k.logp]);
    console.log("sp3: ", data[k.fractioncsp3]);
    // get values from dict: indices 0, 1, 2 correspond to ligands, models, energies
    //let lgdCount = data[0];
    let mdlCount = data[k.modelCounts];
    let energies = data[k.energies];
    let lgdIDs = data[k.ligandIDs];
    let lgdNames = data[k.zincNames];
    let clusters = data[k.clusterAssignment];
    let lgdEnergies = {};
    let lgdNameID = {};

    // build name/ID dict, needed for labelling the small barcharts
    for (let i = 0; i < lgdNames.length; i++) {
      lgdNameID[lgdIDs[i]] = lgdNames[i];
    }

    // fill lgdEnergy dict --> proper format for smallMultiples()
    let nrgIndexCount = 0; // counter to get energy values
    for (let id in lgdIDs) {
      //assumes ligand id always starts with 0
      lgdEnergies[id] = {};
      for (let mdl = 0; mdl < mdlCount[id]; mdl++) {
        // assumes models always start with 0
        lgdEnergies[id][mdl] = {};
        lgdEnergies[id][mdl] = [
          energies[nrgIndexCount],
          mdl,
          clusters[nrgIndexCount],
        ];
        nrgIndexCount++;
      }
    }

    let minVal = d3.min(energies);
    let maxVal = d3.max(energies);

    // gridline function
    let make_gridlines = function (maxY) {
      return d3.axisLeft(yScale).ticks(maxY / 10); //every 10th line
    };

    // create svg object and set width and height attributes
    let svg = d3
      .select("#ligand")
      .append("svg")
      .attr("width", svgWidth + margin.left + margin.right)
      .attr("height", svgHeight + margin.top + margin.bottom - 20)
      .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    /*
		let title = "Ligand binding-energies/AutoDock-scores"
		svg
			.append("text")
			.attr("class", "plotTitle")
			.attr("x", svgWidth / 2 - margin.right/2)
			.attr("y", 0 - margin.top / 2)
			.attr("text-anchor", "middle")
			.text(title);
*/

    // x-axis
    let xScale = d3
      .scaleLinear()
      .domain([Math.ceil(minVal * -1) * -1, Math.floor(maxVal * -1) * -1]) //fixes the smaller bin size of the last bin
      .range([0, svgWidth - margin.left]);
    svg
      .append("g")
      .attr(
        "transform",
        "translate(0," + (svgHeight - margin.bottom - margin.top) + ")"
      )
      .call(d3.axisBottom(xScale));

    // x-label
    svg
      .append("text")
      .attr("class", "xlabel")
      .attr("text-anchor", "start")
      .attr("x", 0)
      .attr("y", svgHeight)
      .style("font-family", "Arial")
      .text("Score / Energy");

    // parameters for the histogram
    let histogram = d3
      .histogram()
      .value(function (d) {
        return d;
      })
      .domain(xScale.domain())
      .thresholds(xScale.ticks(binNum));

    // get bins
    let bins = histogram(energies);

    // y-axis
    let maxValY = d3.max(bins, function (d) {
      return d.length;
    });

    let yScale = d3
      .scaleLinear()
      .domain([0, Math.ceil(maxValY / 10) * 10])
      .range([svgHeight - margin.top - margin.bottom, 0]);

    svg
      .append("g")
      .attr("transform", "translate(" + (svgWidth - margin.left) + ", 0)")
      .call(d3.axisRight(yScale));

    svg
      .append("g")
      .attr("class", "grid")
      .attr("transform", "translate(" + (svgWidth - margin.left) + ", 0)")
      .style("color", "grey")
      .call(
        make_gridlines(maxValY)
          .tickSize(svgWidth - margin.left)
          .tickFormat("")
      );

    // y-label
    svg
      .append("text")
      .attr("class", "ylabel")
      .attr("text-anchor", "end")
      .attr("transform", "rotate(-90)")
      .attr("y", svgWidth)
      .attr("x", 0)
      .style("font-family", "Arial")
      .text("Ligand-Poses");

    // Tooltip

    let tooltip = d3
      .select("#tooltipLigand")
      .style("position", "absolute")
      .style("visibility", "hidden")
      .style("left", margin.left + 10 + "px");

    let mouseOver = function (event, d) {
      d3.select(this).style("stroke", "black").style("stroke-width", "4px");

      tooltip
        .style("visibility", "visible")
        .html(
          "<strong>Range:</strong> " +
            d.x0.toFixed(1) +
            " to " +
            d.x1.toFixed(1) +
            "<br> # " +
            d.length
        )
        .style("left", event.pageX + 10 + "px")
        .style("top", event.pageY - 100 + "px");
    };

    let mouseOut = function () {
      d3.select(this)
        .style("fill", "royalblue")
        .style("opacity", 0.7)
        .style("stroke", "none");
      tooltip.style("visibility", "hidden");
    };

    // add bars
    svg
      .selectAll("rect")
      .data(bins)
      .enter()
      .append("rect")
      .attr("transform", function (d) {
        return "translate(" + xScale(d.x0) + "," + yScale(d.length) + ")";
      }) // shift(left,top) => (x "px", y "px")
      .attr("width", function (d) {
        // d.x0 = lower bound; d.x1 = upper bound
        return Math.abs(xScale(d.x1) - xScale(d.x0) - 1);
      })
      .attr("height", function (d) {
        return svgHeight - yScale(d.length) - margin.bottom - margin.top;
      })
      .style("fill", "royalblue")
      .style("opacity", 0.7)
      .on("mouseover", mouseOver)
      .on("mouseout", mouseOut);

    //smallMultiples(lgdEnergies, minVal, maxVal, lgdNameID, ligandModel);
    return {
      lgdEnergies: lgdEnergies,
      minVal: minVal,
      maxVal: maxVal,
      lgdNameID: lgdNameID,
      ligandModel: ligandModel,
    };
  }
}

/***
 * plots small multiples for all individual ligands
 * @param lgdNameID
 * @param lgdEnergies
 */
export function smallMultiples(plottingOrder, valObj) {
  let margin = { top: 50, right: 10, bottom: 0, left: 40 };

  let box = document.querySelector("#smallMultiples_0").getBoundingClientRect();
  let svgWidth = box.width - (margin.left + margin.right);
  let svgHeight = 200;
  const lgdEnergies = valObj.lgdEnergies;
  const minVal = valObj.minVal;
  const maxVal = valObj.maxVal;
  const lgdNameID = valObj.lgdNameID;

  // create small barchart for each ligand
  for (let i = 0; i < plottingOrder.length; i++) {
    let key = plottingOrder[i];

    let svg = d3
      .select(`#smallMultiples_${i}`)
      .append("svg")
      .attr("width", svgWidth + margin.left + margin.right)
      .attr("height", svgHeight + margin.top + margin.bottom)
      .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    let eList = [];
    let tmpDataList = []; // stores [energy, modelID, assigned Cluster] --> needed to bind model & cluster info to onclick function
    for (let item in lgdEnergies[key]) {
      eList.push(lgdEnergies[key][item][0]);
      tmpDataList.push(lgdEnergies[key][item]);
    }

    // median and mean
    let medianE = d3.median(eList).toFixed(2);
    let meanE = d3.mean(eList).toFixed(2);

    // x-axis
    let modelLabelArray = [];
    Object.keys(lgdEnergies[key]).forEach((element, index) => {
      modelLabelArray.push((index + 1).toString());
    });
    let xScale = d3
      .scaleBand()
      .range([0, svgWidth - margin.left])
      .domain(modelLabelArray)
      .padding(0.4);

    svg.append("g").call(d3.axisTop(xScale));

    // y-axis
    let yScale = d3
      .scaleLinear()
      .domain([Math.ceil(minVal * -1) * -1, Math.floor(maxVal * -1) * -1])
      .range([svgHeight - margin.top - margin.bottom, 0]);
    svg.append("g").call(d3.axisLeft(yScale));

    // gridline function
    let make_gridlines = function (maxY) {
      return d3.axisRight(yScale).ticks(maxY / 0.5); //every 2nd line
    };

    svg
      .append("g")
      .attr("class", "grid")
      .style("color", "grey")
      .call(
        make_gridlines(Math.abs(minVal))
          .tickSize(svgWidth - margin.left)
          .tickFormat("")
      );

    // add title (ligand name)
    svg
      .append("text")
      .attr("x", (svgWidth - margin.left) / 2)
      .attr("y", 0 - margin.top / 2)
      .attr("text-anchor", "middle")
      .style("font-family", "arial")
      .style("font-size", "12px")
      .text(lgdNameID[key] + " (ID: " + key.toString() + ")");

    // add median and mean
    svg
      .append("text")
      .attr("x", (svgWidth - margin.left) / 2)
      .attr("y", svgHeight - margin.top / 2)
      .attr("text-anchor", "middle")
      .style("font-family", "arial")
      .style("font-size", "12px")
      .text("mean: " + meanE + "; " + "median: " + medianE);

    let mouseOver = function () {
      d3.select(this).style("stroke", "black").style("stroke-width", "1px");
    };

    let mouseOut = function () {
      d3.select(this).style("fill", "royalblue").style("stroke", "none");
    };
    // plot bars
    svg
      .selectAll("bars")
      .data(tmpDataList)
      .enter()
      .append("rect")
      .attr("transform", function (d, i) {
        return "translate(" + xScale((i + 1).toString()) + "," + 0 + ")";
      })
      .attr("x", xScale(Object.keys(lgdEnergies[key])))
      .attr("y", yScale(eList))
      .attr("width", xScale.bandwidth)
      .attr("height", function (d) {
        return +yScale(d[0]);
      })
      .style("fill", "royalblue")
      .on("mouseover", mouseOver)
      .on("mouseout", mouseOut)
      .on("click", function (d, i) {
        let mdlID = Number(i[1]);
        let ligandID = Number(key);
        let clusterID = Number(i[2]);
        store.commit("setClusterID", clusterID);
        store.commit("setLigandModelID", [ligandID, mdlID]);

        //store.commit("setSubmitVal", key.toString() +"," + i.toString())
        //document.getElementById("submit").click();
      });
  }
}

/***
 * calculate mean of all ligands and determine order for plotting
 * @type {{}}
 */
export function plot_order(lgdEnergies) {
  let mean_keys = {};
  let plottingOrder = [];

  // create dict with all means
  for (let key in lgdEnergies) {
    let nrg = [];
    for (let item in lgdEnergies[key]) {
      nrg.push(lgdEnergies[key][item][0]);
    }

    mean_keys[key] = d3.mean(nrg).toFixed(2);
  }

  // Create array of keys sorted by mean values (median could be used alternatively)
  // from: https://stackoverflow.com/questions/25500316/sort-a-dictionary-by-value-in-javascript
  let sorted_mean_keys = Object.keys(mean_keys).map(function (key) {
    return [key, mean_keys[key]];
  });

  // Sort the array by medians
  sorted_mean_keys.sort(function (first, second) {
    return first[1] - second[1];
  });

  for (let i = 0; i < sorted_mean_keys.length; i++) {
    plottingOrder.push(parseInt(sorted_mean_keys[i][0]));
  }

  return plottingOrder;
}
