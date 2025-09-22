/*

treat every navigation as a page change that rerenders everything
generate main page
generate library page from a template on nav click to show chart animation

track current lib globally
hold charts globally to destroy on page change

*/

// var currentLib = null;
var currentLib = "RV95";

var charts = [];

google.charts.load("current", { packages: ["corechart", "bar"] });
google.charts.setOnLoadCallback(() => {
  showPage(currentLib);
  window.onresize = () => showPage(currentLib);
});

var pie_chart_default_options = {
  // titleTextStyle: { fontSize: 18 },
  // pieSliceText: "label",
  legend: { position: "left", alignment: "center" },
  // chartArea: { width: "100%", height: "100%" },
  sliceVisibilityThreshold: 0.005, // 0.5%
  pieResidueSliceLabel: "Other",
  animation: { duration: 300, easing: "out" },
};

var main_data = {
  batch_name: "test_batch",
  process_start_time: "2025-07-10T15:18:45.424586-04:00",
  process_end_time: "2025-07-10T15:18:45.692173-04:00",
  current_version: "0.1.0",
  viral_seq_version: "2.0.0",
  number_of_libraries: 2,
  total_reads: 20000,
  raw_sequence_data: [
    ["RV95", 20000],
    ["Other", 2000],
  ],
};

var colors = [
  "#332288",
  "#117733",
  "#44AA99",
  "#88CCEE",
  "#DDCC77",
  "#CC6677",
  "#AA4499",
  "#882255",
];
var bool_colors = { true: colors[3], false: colors[5] };

var lib_data = {
  RV95: {
    raw_distribution: {
      data: [
        ["PR", 5000],
        ["V1V3", 3000],
        ["Other", 2000],
      ],
    },
    raw_sequence_analysis: {
      data: [
        ["GeneralFilterFailed", 2256],
        ["NoMatch", 270],
      ],
      drilldowns: [
        {
          label: "GeneralFilterFailed",
          data: [
            ["IN", 143],
            ["PR", 88],
          ],
        },
        {
          label: "NoMatch",
          data: [
            ["IN", 100],
            ["PR", 100],
          ],
        },
      ],
    },
    number_at_regions: {
      data: [
        ["VIF", 1730, 1711, 1690],
        ["V1V3", 1500, 1600, 1700],
      ],
    },
    detection_sensitivity: {
      data: [
        ["V1V3", 0.3],
        ["PR", 0.5],
      ],
    },
    distinct_to_raw: {
      data: [
        ["IN", 0.007],
        ["PR", 0.003],
      ],
    },
    resampling_index: {
      data: [
        ["IN", 0.007],
        ["PR", 0.003],
      ],
    },
    size_distribution: {
      data: {
        PR: [
          [0, null, 0],
          [0, null, 100],
          [1, 5, null],
          [2, 10, null],
        ],
        IN: [
          [0, null, 0],
          [0, null, 100],
          [1, 5, null],
          [2, 10, null],
        ],
      },
    },
  },
};

// MARK: Navigation
function setNavigation() {
  document.getElementById("nav").innerHTML =
    `
        <a href="https://primer-id.org" target="_BLANK">
            <h3>TCS Log</h3>
        </a>
        <span class="${currentLib == null ? "current_page" : ""}" onclick="showPage(null)">Basic Statistics</span>
    ` +
    Object.keys(lib_data)
      .map(
        (lib, index) =>
          `<span class="${lib === currentLib ? "current_page" : ""}" onclick="showPage('${lib}')">${lib}</span>`,
      )
      .join("");
}

function destroyCharts() {
  charts.forEach((chart) => {
    chart.clearChart();
    google.visualization.events.removeAllListeners(chart);
  });
  charts = [];
}

// MARK: Main Page
function generateMainPage() {
  var pageContent = document.getElementById("page_content");
  var pageBasicStatistics = document.querySelector("#templates .page_main");
  pageContent.innerHTML = pageBasicStatistics.outerHTML;

  // set table fields
  pageContent.querySelector("#batch_name").innerText = main_data["batch_name"];
  pageContent.querySelector("#processed_time").innerText =
    main_data["process_end_time"];
  pageContent.querySelector("#tcs_version").innerText =
    main_data["current_version"];
  pageContent.querySelector("#viral_seq_version").innerText =
    main_data["viral_seq_version"];
  pageContent.querySelector("#number_of_libraries").innerText =
    main_data["number_of_libraries"];
  pageContent.querySelector("#total_reads").innerText =
    main_data["total_reads"];

  // generate chart
  var element = pageContent.querySelector("#raw-sequence-chart");
  var home_chart = new google.visualization.ColumnChart(element);
  home_chart.draw(
    google.visualization.arrayToDataTable([
      ["Library", "Raw Sequences"],
      ...main_data["raw_sequence_data"],
    ]),
    {
      annotations: { alwaysOutside: true },
      //   colors: [bool_colors[true]],
    },
  );
  charts.push(home_chart);
}

// MARK: Lib Page
function generateLibPage(pageID) {
  //set each chart
  var pageContent = document.getElementById("page_content");
  var data = lib_data[pageID];

  pageContent.innerHTML = document.querySelector(
    "#templates .page_lib",
  ).outerHTML;

  // MARK: Raw Distribution
  var raw_distribution = new google.visualization.PieChart(
    pageContent.querySelector("#raw_distribution"),
  );

  let raw_distribution_options = pie_chart_default_options;

  // if (true) {
  //     raw_distribution_options.pieResidueSliceLabel = ""
  // }

  raw_distribution.draw(
    google.visualization.arrayToDataTable([
      ["Region", "Paired Raw"],
      ...data["raw_distribution"]["data"],
    ]),
    raw_distribution_options,
  );
  charts.push(raw_distribution);

  // MARK: Raw Sequence Analysis
  var rsaData = google.visualization.arrayToDataTable([
    ["Categories", "Reason"],
    ...data["raw_sequence_analysis"]["data"],
  ]);
  var raw_sequence_analysis = new google.visualization.PieChart(
    pageContent.querySelector("#raw_sequence_analysis"),
  );
  raw_sequence_analysis.draw(rsaData, pie_chart_default_options);
  charts.push(raw_sequence_analysis);

  // MARK: Drilldowns
  if (data["raw_sequence_analysis"]["drilldowns"].length) {
    var raw_sequence_analysis_drilldown = new google.visualization.PieChart(
      pageContent.querySelector("#raw_sequence_analysis_drilldown"),
    );

    var drilldown_data = data["raw_sequence_analysis"]["drilldowns"][0];
    setDrilldownPieChart(drilldown_data.data, drilldown_data.label);
    charts.push(raw_sequence_analysis_drilldown);

    pageContent.querySelector("#drilldown_label").innerHTML =
      "Raw Sequence Analysis > " + drilldown_data.label;

    // add click event listener to get the key of raw_sequence_analysis pie slice clicked
    google.visualization.events.addListener(
      raw_sequence_analysis,
      "select",
      function () {
        var sel = raw_sequence_analysis.getSelection();

        if (!sel.length || sel[0].row == null) return;

        const r = sel[0].row;
        const label = rsaData.getValue(r, 0);

        const newData = data["raw_sequence_analysis"]["drilldowns"].find(
          (f) => f.label === label,
        );

        setDrilldownPieChart(newData.data, newData.label);
      },
    );

    function setDrilldownPieChart(data, label) {
      pageContent.querySelector("#drilldown_label").innerHTML =
        "Raw Sequence Analysis > " + label;

      raw_sequence_analysis_drilldown.draw(
        google.visualization.arrayToDataTable([
          ["Categories", "Reason"],
          ...data,
        ]),
        {
          title: label,
          ...pie_chart_default_options,
        },
      );
    }
  } else {
    pageContent.querySelector("#drilldown_container").outerHTML = "";
  }

  // MARK: Number at TCS Regions
  var number_at_regions = new google.charts.Bar(
    pageContent.querySelector("#number_at_regions"),
  );
  number_at_regions.draw(
    google.visualization.arrayToDataTable([
      ["Region", "TCS", "Combined TCS", "TCS After QC"],
      ...data["number_at_regions"]["data"],
    ]),
  );
  charts.push(number_at_regions);

  // MARK: Detection Sensitivity
  var detection_sensitivity = new google.visualization.ColumnChart(
    pageContent.querySelector("#detection_sensitivity"),
  );
  detection_sensitivity.draw(
    google.visualization.arrayToDataTable([
      ["Region", "Detection Sensitivity"],
      ...data["detection_sensitivity"]["data"],
    ]),
    {
      legend: {
        position: "none",
      },
      vAxis: {
        scaleType: "log",
        viewWindow: {
          min: 0.00001,
        },
        ticks: [0.00001, 0.0001, 0.001, 0.01, 0.1, 1],
      },
    },
  );
  charts.push(detection_sensitivity);

  // MARK: Distinct To Raw
  var distinct_to_raw = new google.visualization.ColumnChart(
    pageContent.querySelector("#distinct_to_raw"),
  );
  distinct_to_raw.draw(
    google.visualization.arrayToDataTable([
      ["Region", "Distinct to Raw"],
      ...data["distinct_to_raw"]["data"],
    ]),
    {
      legend: {
        position: "none",
      },
    },
  );
  charts.push(distinct_to_raw);

  // MARK: Resampling Index
  var resampling_index = new google.visualization.ColumnChart(
    pageContent.querySelector("#resampling_index"),
  );
  resampling_index.draw(
    google.visualization.arrayToDataTable([
      ["Region", "Resampling Index"],
      ...data["resampling_index"]["data"],
    ]),
    {
      legend: {
        position: "none",
      },
    },
  );
  charts.push(resampling_index);

  // MARK: Size Distribution
  var size_distribution_container =
    pageContent.querySelector("#size_distribution");

  Object.keys(data["size_distribution"]["data"]).map((region) => {
    var regionData = data["size_distribution"]["data"][region];

    var region_container = document.createElement("div");
    size_distribution_container.appendChild(region_container);

    var size_distribution = new google.visualization.ComboChart(
      region_container,
    );

    var view = new google.visualization.DataView(
      google.visualization.arrayToDataTable([
        ["Index", "Distribution", "Cutoff"],
        ...regionData,
      ]),
    );
    size_distribution.draw(view, {
      pointSize: 5,
      title: region,
      hAxis: { title: "Raw sequencing reads per unique PID" },
      vAxis: {
        title: "# of PIDs ",
        logScale: true,
      },
      // colors: ["' + bool_colors[true] + '"],
      legend: "none",
      seriesType: "scatter",
      series: {
        1: {
          type: "line",
          // color: "' + bool_colors[false] + '",
          pointsVisible: false,
        },
      },
    });

    charts.push(size_distribution);
  });
}

function showPage(pageID) {
  destroyCharts();

  document.getElementById("chart_error_message").innerHTML = "";

  currentLib = pageID;

  console.log({ pageID });

  try {
    if (pageID == null) {
      generateMainPage();
    } else {
      generateLibPage(pageID);
    }
  } catch (e) {
    console.log(e);
    document.getElementById("chart_error_message").innerHTML =
      "There was an error displaying your data.<br>To help us improve, please forward this file or link to<br>clarkmu@email.unc.edu";
  }

  setNavigation(pageID);
}

// MARK: Coloring

var regionColors = {
  PR: "#1f77b4", // blue
  IN: "#ff7f0e", // orange
  V1V3: "#2ca02c", // green
  VIF: "#9467bd", // purple
  // "Other" not fixed → falls back to random
};

function getRegionColor(region) {
  if (regionColors[region]) return regionColors[region];
  if (region === "Other") {
    // cache one random value so "Other" looks the same across charts
    if (!regionColors.Other) regionColors.Other = randomHex();
    return regionColors.Other;
  }
  // unknown region fallback
  if (!regionColors[region]) regionColors[region] = randomHex();
  return regionColors[region];
}

function hexToHsl(hex) {
  hex = hex.replace(/^#/, "");
  if (hex.length === 3) {
    hex = hex
      .split("")
      .map((ch) => ch + ch)
      .join("");
  }
  const r = parseInt(hex.substring(0, 2), 16) / 255;
  const g = parseInt(hex.substring(2, 4), 16) / 255;
  const b = parseInt(hex.substring(4, 6), 16) / 255;

  const max = Math.max(r, g, b),
    min = Math.min(r, g, b);
  let h,
    s,
    l = (max + min) / 2;

  if (max === min) {
    h = s = 0; // achromatic
  } else {
    const d = max - min;
    s = l > 0.5 ? d / (2 - max - min) : d / (max + min);
    switch (max) {
      case r:
        h = (g - b) / d + (g < b ? 6 : 0);
        break;
      case g:
        h = (b - r) / d + 2;
        break;
      case b:
        h = (r - g) / d + 4;
        break;
    }
    h /= 6;
  }
  return { h: h * 360, s: s, l: l };
}

function gradientShades(hex, steps = 4) {
  const hsl = hexToHsl(hex);
  const shades = [];
  for (let i = 0; i < steps; i++) {
    const factor = (i / (steps - 1)) * 0.5; // 0 → base, 0.5 → lighter
    shades.push(d3.hsl(hsl.h, hsl.s, Math.min(1, hsl.l + factor)).formatHex());
  }
  return shades;
}
