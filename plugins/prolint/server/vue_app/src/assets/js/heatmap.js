import * as d3 from "d3";
import { clusterDataKeys as k, store } from "../../main.js";

export function drawColormap(in_margin, svgWidth, values, state) {
  let margin = { top: 0, right: 0, bottom: 0, left: 0 };
  margin = in_margin;
  margin.top = 20;
  svgWidth = svgWidth - 20;

  let bars = d3
    .select("#colormap")
    .append("svg")
    .attr("width", svgWidth / 2 + margin.left + margin.right)
    .attr("height", 20 + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  // inspired by: https://stackoverflow.com/questions/20847161/how-can-i-generate-as-many-colors-as-i-want-using-d3/30912617

  const svgHeight = 25;

  // Append title
  let title = "Score / Energy";
  //let title = "Binding-Energy / AutoDock-Score";
  bars
    .append("text")
    .attr("transform", "translate(" + (margin.left + -77) + ",-7)")
    .attr("class", "plotTitle")
    .attr("x", margin.left + margin.right)
    .attr("y", 0)
    .attr("text-anchor", "middle")
    .text(title);

  bars
    .selectAll()
    .data(d3.range(0, 1, 1 / state.n))
    .enter()
    .append("rect")
    .attr("class", "bars")
    .attr("y", 0)
    .attr("height", svgHeight / 2)
    .attr("width", svgWidth / 2 / state.n + 1)
    .attr("x", (d, i) => (i * (svgWidth / 2)) / state.n)
    .style("fill", state.colorScale);

  let x = d3
    .scaleLinear()
    .range([0, (svgWidth / 2 / state.n) * state.n])
    .domain([
      store.state.dynamicEnergyMin ? store.state.minEnergy : 0,
      store.state.maxEnergy,
    ]);

  bars
    .append("g")
    .attr("class", "axis")
    .attr("transform", "translate(0," + (svgHeight / 2 + 5) + ")")
    .call(
      d3
        .axisBottom(x)
        .ticks(6)
        .tickSize(5)
        .tickFormat((d) => "-" + d)
    )
    .select(".domain");

  // push empty object to initiate recalculation of SVGs
  store.commit("setColorMapSVGs", {});
  //bars.exit().remove();
}

export function drawColormapTooltip(in_margin, svgHeight, val) {
  let margin = { top: 0, right: 0, bottom: 0, left: 0 };
  margin = in_margin;

  const state = {
    rangeStart: 0,
    n: 200,
    colorScale: d3.scaleSequential(d3.interpolateInferno),
  };

  //let bars = d3.select("#colormap_tooltip");
  //bars.select("svg").empty();

  //d3.select("#colormap_tooltip").select("svg").selectAll("*").remove();
  //d3.selectAll("#ccc").remove();

  //console.log("poss:", XX, XX.offsetParent.offsetLeft+XX.offsetLeft);
  //let xX = XX.offsetTop;
  //let yY = XX.offsetParent.offsetLeft+XX.offsetLeft+100;
  const svgWidth = 30;

  let barsSVG = d3
    .create("svg")
    .style("z-index", "200")
    .attr("width", svgWidth / 2 + margin.left + margin.right)
    .attr("height", svgHeight / 2 + margin.top + margin.bottom);

  let bars = barsSVG
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  // inspired by: https://stackoverflow.com/questions/20847161/how-can-i-generate-as-many-colors-as-i-want-using-d3/30912617

  // Append title
  let title = "Score / Energy";
  //let title = "Binding-Energy / AutoDock-Score";
  bars
    .append("text")
    .attr("transform", "rotate(-90), translate(" + (0 - 225) + ",-40)")
    .attr("class", "plotTitle")
    .attr("x", 0)
    .attr("y", 0)
    .attr("text-anchor", "end")
    .style("font-size", "12px")
    .text(title);

  bars
    .selectAll()
    .data(d3.range(1, 0, -1 / state.n))
    .enter()
    .append("rect")
    .attr("class", "bars")
    .attr("y", (d, i) => (i * (svgHeight / 2)) / state.n)
    .attr("height", svgHeight / 2 / state.n + 1)
    .attr("width", svgWidth / 2)
    .attr("x", 0)
    .style("fill", state.colorScale);

  let x = d3
    .scaleLinear()
    .range([(svgHeight / 2 / state.n) * state.n, 0])
    .domain([
      store.state.dynamicEnergyMin ? store.state.minEnergy : 0,
      store.state.maxEnergy,
    ]);

  bars
    .append("g")
    .attr("class", "axis")
    .attr("transform", "translate(" + -5 + "," + 0 + ")")
    .call(
      d3
        .axisLeft(x)
        .ticks(6)
        .tickFormat((d) => "-" + d)
        .tickSize(5)
    )
    .style("color", "#000000");
  //.select(".domain");

  var color = "green";
  var triangleSize = 100;
  //var verticalTransform = 10 + Math.sqrt(triangleSize);

  var triangle = d3.symbol().type(d3.symbolTriangle).size(triangleSize);
  let scaled = x(Number(val * -1));

  bars
    .append("g")
    .append("path")
    .attr("d", triangle)
    .attr("stroke", color)
    .attr("fill", color)
    .attr("transform", "translate(" + 20 + "," + scaled + ")rotate(-90)");

  //bars.exit().remove();
  var serializer = new XMLSerializer();
  var source = serializer.serializeToString(barsSVG.node());

  return source;
}

export function drawHeatmap(dataset, div, axisDiv) {
  d3.select(div).selectAll("svg").remove();
  d3.select(axisDiv).selectAll("svg").remove();
  d3.select("#colormap").selectAll("svg").remove();

  // for drawColormap()
  const colorMapState = {
    rangeStart: 0,
    n: 200,
    colorScale: d3.scaleSequential(d3.interpolateInferno),
  };

  // Check if data available
  if (Object.keys(dataset).length === 0 && dataset.constructor === Object) {
    console.log("No data available. Please send data first.");
  } else {
    let withNoise = store.state.showNoise;
    // inspired by: https://www.d3-graph-gallery.com/graph/heatmap_style.html
    let clientWidth = document.getElementById("heatmap").clientWidth;
    //let clientHeight = document.getElementById("heatmap").clientHeight;

    // set clusterCnt (with/without noise)
    let clustIDCnt;
    if (withNoise) {
      clustIDCnt = dataset[k.clusterSizes].length + 1;
    } else {
      clustIDCnt = dataset[k.clusterSizes].length;
    }
    console.log("clustIDCnt", clustIDCnt);

    var margin = { top: 0, right: 30, bottom: 30, left: 50 };

    let gap = 1;
    const fontSizePx =
      parseFloat(getComputedStyle(document.documentElement).fontSize) * 2;
    clientWidth = clientWidth - margin.left - margin.right - gap * clustIDCnt;
    let rectSize = clientWidth / clustIDCnt;
    let svgWidth = clustIDCnt * (rectSize + gap);
    let svgHeight = dataset[3].length * (fontSizePx + gap);
    let xLabelSpace = 19;

    let addNoise = 0;
    if (withNoise) {
      addNoise = 1;
    }
    console.log("addNoise", addNoise);

    d3.select("#vueSwitch_showNoise")
      .style("margin-left", margin.left + "px")
      .style("margin-bottom", -margin.top + 10 + "px");

    let axisSvg = d3
      .select("#heatmapAxis")
      .append("svg")
      .attr("width", svgWidth + margin.left + margin.right)
      .attr("height", 30)
      .append("g")
      .attr("transform", "translate(" + margin.left + "," + 0 + ")");

    // Append the svg object
    let svg = d3
      .select(div)
      .append("svg")
      .attr("width", svgWidth + margin.left + margin.right)
      .attr("height", svgHeight + margin.top + margin.bottom)
      .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    // Labels of row and columns -> unique identifier of the column called 'group' and 'variable'
    let mdlsPerLigand,
      energies = [],
      ligandIDs,
      gblMdlID_clusterAssignment = [];
    mdlsPerLigand = dataset[k.modelCounts];
    energies = dataset[k.energies];
    ligandIDs = dataset[k.ligandIDs];
    gblMdlID_clusterAssignment = dataset[k.clusterAssignment];

    drawColormap(margin, svgWidth, energies, colorMapState);

    let totalSubSquareCounter = [],
      emptySquares = [];
    let subRectCounter = [];
    let cnt = 0;

    // initialize array with true
    for (let i = 0; i < ligandIDs.length; i++) {
      for (let j = 0; j < clustIDCnt; j++) {
        totalSubSquareCounter.push([]);
        totalSubSquareCounter[i][j] = 0;
        subRectCounter.push([]);
        subRectCounter[i][j] = 0;
      }
    }

    cnt = 0;

    for (let i = 0; i < ligandIDs.length; i++) {
      for (let j = 0; j < mdlsPerLigand[i]; j++) {
        if (gblMdlID_clusterAssignment[cnt] != -1 || withNoise) {
          totalSubSquareCounter[i][
            gblMdlID_clusterAssignment[cnt] + addNoise
          ]++;
          //console.log(i,j,gblMdlID_clusterAssignment[ligandModelIDs.length-1]+addNoise)
        }
        cnt++;
      }
    }

    // determine sub rects
    cnt = 0;
    let all = [];
    for (let i = 0; i < ligandIDs.length; i++) {
      for (let j = 0; j < mdlsPerLigand[i]; j++) {
        if (gblMdlID_clusterAssignment[cnt] != -1 || withNoise) {
          let struct = {
            ligandID: 0,
            modelID: 0,
            gblMdlID: 0,
            clusterID: 0,
            internIndex: 0,
            totalCnt: 0,
          };
          struct.ligandID = i;
          struct.modelID = j;
          struct.gblMdlID = cnt;
          struct.clusterID = gblMdlID_clusterAssignment[cnt];
          struct.internIndex =
            subRectCounter[i][gblMdlID_clusterAssignment[cnt] + addNoise];
          //console.log("intern: ",subRectCounter[i][gblMdlID_clusterAssignment[cnt]+addNoise],"total: ",totalSubSquareCounter[i][gblMdlID_clusterAssignment[cnt]+addNoise])
          subRectCounter[i][gblMdlID_clusterAssignment[cnt] + addNoise]++;
          struct.totalCnt =
            totalSubSquareCounter[i][
              gblMdlID_clusterAssignment[cnt] + addNoise
            ];
          all.push(struct);
          //console.log(i,j,gblMdlID_clusterAssignment[ligandModelIDs.length-1]+addNoise)
        }
        cnt++;
      }
    }

    cnt = 0;
    let start = 0;
    if (withNoise) {
      start = -1;
    }
    for (let i = 0; i < ligandIDs.length; i++) {
      for (let j = start; j < clustIDCnt; j++) {
        let emptyArray = [];
        if (totalSubSquareCounter[i][j] == 0) {
          emptySquares.push(emptyArray);
          emptySquares[cnt].push(i, j);
          cnt++;
        }
      }
    }

    let xDomain = [
      0,
      d3.max(
        gblMdlID_clusterAssignment.map(function (d) {
          return d;
        })
      ),
    ];
    if (withNoise) {
      xDomain = [-1, xDomain[1]];
    }

    // Build X scales and axis:
    let x = d3
      .scaleBand()
      .range([0, svgWidth])
      .domain(d3.range(xDomain[0], xDomain[1] + 1, 1))
      .padding(0.0);
    axisSvg
      .append("g")
      .attr("class", "axis")
      .attr("transform", "translate(0," + 15 + ")")
      .call(
        d3
          .axisBottom(x)
          .tickSize(0)
          .tickFormat(function (d) {
            if (d != -1) {
              return d;
            } else {
              return "noise";
            }
          })
      )
      .select(".domain")
      .remove();

    svg
      .append("g")
      .attr("class", "axis")
      .attr("transform", "translate(0," + svgHeight + ")")
      .call(
        d3
          .axisBottom(x)
          .tickSize(0)
          .tickFormat(function (d) {
            if (d != -1) {
              return d;
            } else {
              return "noise";
            }
          })
      )
      .select(".domain")
      .remove();

    // add x axis text-label TOP
    axisSvg
      .append("text")
      .attr("transform", "translate(0," + 10 + ")")
      .attr("text-anchor", "left")
      .style("font-size", "14px")
      .style("fill", "black")
      .text("Cluster ID");

    // make positioning of y lable depend on length of ligandID ("128" = length of 3)
    let tmp = ligandIDs[ligandIDs.length - 1].toString().length;
    console.log("IDLength: ", tmp);
    let yLabelSpace_top = tmp * 10 + 3;
    tmp = ligandIDs[0].toString().length;
    let yLabelSpace_bot = tmp * 10 + 5;
    // add y axis text-label TOP
    svg
      .append("text")
      .attr("transform", "rotate(-90),translate(-65,-" + yLabelSpace_top + ")")
      .attr("text-anchor", "left")
      .style("font-size", "14px")
      .style("fill", "black")
      .text("Ligand ID");

    // add x axis text-label BOTTOM
    svg
      .append("text")
      .attr("transform", "translate(0," + (svgHeight + xLabelSpace + 10) + ")")
      .attr("text-anchor", "left")
      .style("font-size", "14px")
      .style("fill", "black")
      .text("Cluster ID");

    // add y axis text-label BOTTOM
    svg
      .append("text")
      .attr(
        "transform",
        "rotate(-90), translate(-" + svgHeight + ",-" + yLabelSpace_bot + ")"
      )
      .attr("text-anchor", "left")
      .style("font-size", "14px")
      .style("fill", "black")
      .text("Ligand ID");

    // Build Y scales and axis:
    let y = d3.scaleBand().range([svgHeight, 0]).domain(ligandIDs).padding(0.0);
    svg
      .append("g")
      .attr("class", "axis")
      .call(d3.axisLeft(y).tickSize(0))
      .select(".domain")
      .remove();

    // create a tooltip
    let tooltip = d3
      .select("#tooltipHeatmap")
      .style("position", "absolute")
      .style("opacity", 0)
      .attr("class", "tooltip");

    // Three function that change the tooltip when user hover / move / leave a cell
    let mouseover = function (event) {
      tooltip.style("opacity", 1);
      d3.select(this)
        .style("stroke-width", 3)
        .style("stroke", "black")
        .style("opacity", 1)
        .style("left", event.pageX + 10 + "px")
        .style("top", event.pageY - 100 + "px");
    };

    let mousemove = function (event, d) {
      //const e = d3.select(this).nodes();
      //const i = e.indexOf(this);

      let gblMdlID = d["gblMdlID"];
      tooltip
        .html(
          "<center><h1>" +
            energies[gblMdlID].toFixed(1) +
            "</h1></center>" +
            "<strong>" +
            dataset[k.zincNames][d["ligandID"]] +
            "</strong>" +
            "<br>" +
            "<table>" +
            "<tr><td> <strong>Cluster:</strong></td><td> " +
            d["clusterID"] +
            "</td></tr>" +
            "<tr><td><strong>Ligand.:</strong></td><td> " +
            d["ligandID"] +
            "</td></tr>" +
            "<tr><td><strong>Model.:</strong></td><td> " +
            (parseInt(d["modelID"]) + 1) +
            "</td></tr>" +
            "</table>"
        )
        .style("left", event.pageX + 10 + "px")
        .style("top", event.pageY - 100 + "px");
    };

    let mouseleave = function (event) {
      tooltip
        .style("opacity", 0)
        .style("left", event.pageX + 10 + "px")
        .style("top", event.pageY - 100 + "px");
      d3.select(this)
        .style("stroke-width", 0)
        .style("stroke", "white")
        .style("opacity", 0.8);
    };

    let mouseclick = function (event, d) {
      //let gblMdlID = d["gblMdlID"];
      //store.commit("setGblMdlID", gblMdlID);
      console.log('d["ligandID"]', d["ligandID"], 'd["modelID"]', d["modelID"]);
      store.commit("setClusterID", d["clusterID"]);
      //store.commit("setLigandID", [d["ligandID"]]);
      store.commit("setLigandModelID", [d["ligandID"], [d["modelID"]]]);
    };

    let roundness = 2;
    let value = 0;
    // add the squares
    svg
      .selectAll()
      .data(all, function (d, i) {
        return i;
      })
      .enter()
      .append("rect")
      .attr("x", function (d) {
        return (
          x(gblMdlID_clusterAssignment[d["gblMdlID"]]) +
          (rectSize / d["totalCnt"]) * d["internIndex"] +
          gap
        );
      })
      .attr("y", function (d) {
        return y(d["ligandID"]) + gap;
      })
      .attr("rx", roundness)
      .attr("ry", roundness)
      .attr("width", function (d) {
        return Math.floor(rectSize / d["totalCnt"]);
      })
      .attr("height", fontSizePx)
      .style("fill", function (d) {
        return store.getters.colorScale(energies[d["gblMdlID"]]);
      })
      .attr("score", value)
      .style("stroke-width", 0)
      .style("stroke", "white")
      .style("opacity", 0.8)
      .on("mouseover", mouseover)
      .on("mousemove", mousemove)
      .on("mouseleave", mouseleave)
      .on("click", mouseclick);

    // grey rects for empty fields
    svg
      .selectAll()
      .data(emptySquares, function (event, d) {
        return d;
      })
      .enter()
      .append("rect")
      .attr("x", function (event, d) {
        return x(emptySquares[d][1] - addNoise) + gap;
      })
      .attr("y", function (event, d) {
        return y(emptySquares[d][0]) + gap;
      })
      .attr("rx", roundness)
      .attr("ry", roundness)
      .attr("width", rectSize)
      .attr("height", fontSizePx)
      .style("fill", d3.color("rgba(70%,70%,70%,0.8)").formatHex())
      .attr("score", value)
      .style("stroke-width", 4)
      .style("stroke", "none")
      .style("opacity", 0.8);

    // Add subtitle to graph
    svg
      .append("text")
      .attr("x", 0)
      .attr("y", -20)
      .attr("text-anchor", "left")
      .style("font-size", "14px")
      .style("fill", "grey")
      .style("max-width", 400)
      .text("");
  }
}
