import * as d3 from "d3";
import { store } from "../../main.js";

export class functionalGroupChart {
  margin = { top: 20, left: 10, bottom: 10, right: 10 };
  horizontalAxisSvg;
  horizontalAxis;
  horizontalScale;

  mainSvg;
  mainSvgWidth;

  barChart;
  barHeight = 23;
  barPadding = 4;
  amtBars = 0;

  divs; // 0 horizontal axis, 1 mainSVG

  mainData;
  mainDataExtent;

  constructor(divs, data) {
    this.divs = divs;
    //this.mainData = Object.keys(data).map((key) => [Number(key), data[key]]);
    this.mainData = data.map((dat) => [dat.funcGroupID, dat.funcGroupCnt]);
    this.amtBars = data.length;
    this.mainData.sort((a, b) => d3.descending(a[1], b[1]));
    this.init();
    console.log("funcGroupData:", data);
  }

  init() {
    this.mainDataExtent = [
      0,
      Math.max(...this.mainData.map((elem) => elem[1])),
    ];
    this.mainSvgWidth = document.getElementById(
      this.divs[1].replace("#", "")
    ).clientWidth;
    console.log("chart width", this.mainSvgWidth);
    this.horizontalScale = d3
      .scaleLinear()
      .domain(this.mainDataExtent)
      .nice()
      .range([this.margin.left, this.mainSvgWidth - this.margin.right]);
    this.drawBars();
    this.drawHorizontalAxis();
    // Tooltip
    this.tooltip = d3.select("#funcGroupTooltip");
    this.tooltip
      .style("position", "absolute")
      .style("visibility", "hidden")
      .style("left", this.margin.left + 10 + "px");
  }
  updateData(newData) {
    //this.mainData =  Object.keys(newData).map((key) => [Number(key), newData[key]]);
    this.mainData = newData.map((dat) => [dat.funcGroupID, dat.funcGroupCnt]);
    this.amtBars = this.mainData.length;
    //this.mainSvg.attr("height", this.margin.top + this.margin.bottom + this.amtBars * (this.barHeight+this.barPadding))
    this.mainDataExtent = [
      0,
      Math.max(...this.mainData.map((elem) => elem[1])),
    ];
    this.mainData.sort((a, b) => d3.descending(a[1], b[1]));
    this.horizontalScale.domain(this.mainDataExtent);

    this.barDrawLoop(this.mainData);
    this.updateHorizontalAxis();
  }
  drawBars() {
    // inspired by https://observablehq.com/@d3/bar-chart-race and https://observablehq.com/@d3/hierarchical-bar-chart
    let selectedDiv = d3.select(this.divs[1]);
    this.mainSvg = selectedDiv
      .append("svg")
      .attr("id", "mainsvg")
      .attr("width", this.mainSvgWidth)
      .attr(
        "height",
        this.amtBars * (this.barHeight + this.barPadding) +
          this.margin.top +
          this.margin.bottom
      );
    this.mainG = this.mainSvg.append("g");
    this.barDrawLoop(this.mainData);
  }
  barDrawLoop(data) {
    this.mainSvg
      .selectAll("g")
      .data(data, (d) => d)
      .join(
        (enter) => this.barsEnter(enter),
        (update) => update, //this.barsUpdate(update),
        (exit) => this.barsExit(exit)
      );
  }
  barsEnter(selection) {
    let currG = selection.append("g");
    let self = this;
    currG
      .append("rect")
      .attr("x", this.horizontalScale(0))
      .attr("y", (d, i) => i * (this.barHeight + this.barPadding))
      .attr("width", 0)
      .attr("height", this.barHeight)
      .attr("stroke", "black")
      .attr("stroke-width", 1)
      .attr("fill", "royalblue")
      .attr("opacity", 0.5)
      .attr("cursor", "pointer")
      .on("click", (event, d) => store.commit("setSelectedFuncGroup", d[0]))
      .call((enter) =>
        enter
          .transition(this.mainSvg.transition().duration(750))
          .attr(
            "width",
            (d) => this.horizontalScale(d[1]) - this.horizontalScale(0)
          )
      )
      .on("mouseover", function (event, d) {
        const data = d;
        self.tooltip
          .style("visibility", "visible")
          .html("<br>" + data[0] + "<br><strong>" + data[1] + "</strong>")
          .style("left", event.layerX + "px")
          .style("top", event.layerY + "px");
      });

    // for textlabels https://observablehq.com/@d3/horizontal-bar-chart
    currG
      .append("text")
      .attr("x", (d) => this.horizontalScale(d[1]))
      .attr(
        "y",
        (d, i) => i * (this.barHeight + this.barPadding) + this.barHeight / 2
      )
      .attr("dx", -4)
      .attr("dy", "0.35em")
      .text((d) => store.state.functionalGroupsIDsToWord[d[0]])
      .attr("text-anchor", "end")
      .call((text) =>
        text
          .filter(
            (d) =>
              this.horizontalScale(d[1]) +
                String(store.state.functionalGroupsIDsToWord[d[0]]).length *
                  10 <
              this.mainSvgWidth
          ) // short bars
          .attr("dx", +4)
          .attr("fill", "black")
          .attr("text-anchor", "start")
      );
  }
  barsUpdate(selection) {
    let bars = selection.selectAll("rect");
    let text = selection.selectAll("text");
    bars
      .attr("x", this.horizontalScale(0))
      .attr(
        "width",
        (d) => this.horizontalScale(d[1]) - this.horizontalScale(0)
      )
      .attr("height", this.barHeight)
      .attr("stroke", "black")
      .attr("stroke-width", 1)
      .attr("fill", "green")
      .attr("opacity", 0.5);
    text
      .attr("x", (d) => this.horizontalScale(d[1]) - String(d[0]).length * 12)
      .text((d) => store.state.functionalGroupsIDsToWord[d[0]])
      .call((text) =>
        text
          .append("tspan")
          .attr("fill-opacity", 0.7)
          .attr("font-weight", "normal")
      );
    //.attr("x", -6)
    //.attr("dy", "1.15em"))
  }
  barsExit(selection) {
    let bars = selection.selectAll("rect");
    let text = selection.selectAll("text");

    bars
      .attr("fill", "brown")
      .call((exit) =>
        exit
          .transition(this.mainSvg.transition().duration(750))
          .attr("width", 0)
          .remove()
      );
    text
      .attr("fill", "brown")
      .call((exit) =>
        exit
          .transition(this.mainSvg.transition().duration(750))
          .attr("opacity", 0)
          .remove()
      );

    selection.call((exit) =>
      exit.transition(this.mainSvg.transition().duration(750)).remove()
    );
  }
  drawHorizontalAxis() {
    let selectedDiv = d3.select(this.divs[0]);
    this.horizontalAxisSvg = selectedDiv
      .append("svg")
      .attr("width", this.mainSvgWidth)
      .attr("height", this.margin.top + this.margin.bottom);
    let g = this.horizontalAxisSvg
      .append("g")
      .attr("id", "topaxis")
      .attr("transform", `translate(0,${this.margin.top})`);
    this.horizontalAxis = g.call(d3.axisTop(this.horizontalScale));
  }
  updateHorizontalAxis() {
    this.horizontalAxisSvg
      .select("#topaxis")
      .transition()
      .duration(750)
      .call(d3.axisTop(this.horizontalScale));
  }
}
