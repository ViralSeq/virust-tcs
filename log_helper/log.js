/*

treat every navigation as a page change that rerenders everything
generate main page
generate library page from a template on nav click to show chart animation

track current lib globally
hold charts globally to destroy on page change

*/

var currentLib = null;
// var currentLib = "TEST_DATA";

let { main_data, lib_data } = app_data;

// insert test data if no data
if (!main_data.batch_name) {
  main_data = {
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

  lib_data = {
    TEST_DATA: {
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
          ["VIF", 1730, null, 1711],
          ["V1V3", 1500, 1600, 1700],
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
}

var charts = [];

google.charts.load("current", { packages: ["corechart", "bar"] });
google.charts.setOnLoadCallback(() => {
  showPage(currentLib);
  window.onresize = () => showPage(currentLib);
});

var pie_chart_default_options = {
  legend: { position: "left", alignment: "center" },
  sliceVisibilityThreshold: 0.005,
  pieResidueSliceLabel: "Other",
  animation: { duration: 300, easing: "out" },
  pieHole: 0.4,
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

  var element = pageContent.querySelector("#raw-sequence-chart");
  var home_chart = new google.visualization.ColumnChart(element);
  home_chart.draw(
    google.visualization.arrayToDataTable([
      ["Library", "Raw Sequences"],
      ...main_data["raw_sequence_data"],
    ]),
    {
      annotations: { alwaysOutside: true },
    },
  );
  charts.push(home_chart);
}

// MARK: Lib Page
function generateLibPage(pageID) {
  var pageContent = document.getElementById("page_content");
  var data = lib_data[pageID];

  pageContent.innerHTML = document.querySelector(
    "#templates .page_lib",
  ).outerHTML;

  var raw_distribution = new google.visualization.PieChart(
    pageContent.querySelector("#raw_distribution"),
  );

  let raw_distribution_options = {
    ...pie_chart_default_options,
    colors: data["raw_distribution"]["data"].map(function (r) {
      return getRegionColor(r[0]);
    }),
  };

  raw_distribution.draw(
    google.visualization.arrayToDataTable([
      ["Region", "Paired Raw"],
      ...data["raw_distribution"]["data"],
    ]),
    raw_distribution_options,
  );
  charts.push(raw_distribution);

  var rsaData = google.visualization.arrayToDataTable([
    ["Categories", "Reason"],
    ...data["raw_sequence_analysis"]["data"],
  ]);
  var raw_sequence_analysis = new google.visualization.PieChart(
    pageContent.querySelector("#raw_sequence_analysis"),
  );
  var parentLabels = data["raw_sequence_analysis"]["data"].map(function (r) {
    return r[0];
  });
  var parentColorByLabel = {};
  parentLabels.forEach(function (label) {
    parentColorByLabel[label] = getRegionColor(label);
  });
  raw_sequence_analysis.draw(rsaData, {
    ...pie_chart_default_options,
    colors: parentLabels.map(function (l) {
      return parentColorByLabel[l];
    }),
  });
  charts.push(raw_sequence_analysis);

  if (data["raw_sequence_analysis"]["drilldowns"].length) {
    var raw_sequence_analysis_drilldown = new google.visualization.PieChart(
      pageContent.querySelector("#raw_sequence_analysis_drilldown"),
    );

    var drilldown_data = data["raw_sequence_analysis"]["drilldowns"][0];
    setDrilldownPieChart(drilldown_data.data, drilldown_data.label);
    charts.push(raw_sequence_analysis_drilldown);

    pageContent.querySelector("#drilldown_label").innerHTML =
      "Raw Sequence Analysis > " + drilldown_data.label;

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

    function setDrilldownPieChart(dataRows, label) {
      pageContent.querySelector("#drilldown_label").innerHTML =
        "Raw Sequence Analysis > " + label;

      var base = parentColorByLabel[label] || "#8888ee";
      var shades = gradientShades(base, dataRows.length || 2);

      raw_sequence_analysis_drilldown.draw(
        google.visualization.arrayToDataTable([
          ["Categories", "Reason"],
          ...dataRows,
        ]),
        {
          title: label,
          ...pie_chart_default_options,
          colors: shades,
        },
      );
    }
  } else {
    pageContent.querySelector("#drilldown_container").outerHTML = "";
  }

  var number_at_regions = new google.charts.Bar(
    pageContent.querySelector("#number_at_regions"),
  );
  number_at_regions.draw(
    google.visualization.arrayToDataTable([
      ["Region", "TCS", "Combined TCS", "TCS After QC"],
      ...data["number_at_regions"]["data"],
    ]),
    {
      colors: [ggplotPalette[0], ggplotPalette[1], ggplotPalette[2]],
    },
  );
  charts.push(number_at_regions);

  var distinct_to_raw = new google.visualization.ColumnChart(
    pageContent.querySelector("#distinct_to_raw"),
  );
  var dtrRows = [["Region", "Distinct to Raw", { role: "style" }]].concat(
    data["distinct_to_raw"]["data"].map(function (r) {
      return [r[0], r[1], getRegionColor(r[0])];
    }),
  );
  distinct_to_raw.draw(google.visualization.arrayToDataTable(dtrRows), {
    legend: {
      position: "none",
    },
  });
  charts.push(distinct_to_raw);

  var resampling_index = new google.visualization.ColumnChart(
    pageContent.querySelector("#resampling_index"),
  );
  var riRows = [["Region", "Resampling Index", { role: "style" }]].concat(
    data["resampling_index"]["data"].map(function (r) {
      return [r[0], r[1], getRegionColor(r[0])];
    }),
  );
  resampling_index.draw(google.visualization.arrayToDataTable(riRows), {
    legend: {
      position: "none",
    },
  });
  charts.push(resampling_index);

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
    var base = getRegionColor(region);
    var shade = gradientShades(base, 2)[1];
    size_distribution.draw(view, {
      pointSize: 5,
      title: region,
      hAxis: { title: "Raw sequencing reads per unique PID" },
      vAxis: {
        title: "# of PIDs ",
        logScale: true,
      },
      legend: "none",
      seriesType: "scatter",
      series: {
        0: { color: base },
        1: { type: "line", pointsVisible: false, color: shade },
      },
    });

    charts.push(size_distribution);
  });
}

function showPage(pageID) {
  destroyCharts();

  document.getElementById("chart_error_message").innerHTML = "";

  currentLib = pageID;

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
  PR: "#F8766D", // red
  IN: "#7CAE00", // green
  V1V3: "#00BFC4", // cyan/blue
  VIF: "#C77CFF", // purple
  GAG: "#E68613", // orange
  NEF: "#00A9FF", // light blue
  POL: "#52B415", // lime
  ENV: "#FF61CC", // pink
  Other: "#999999", // gray fallback
};

var ggplotPalette = [
  "#F8766D",
  "#7CAE00",
  "#00BFC4",
  "#C77CFF",
  "#E68613",
  "#00A9FF",
  "#52B415",
  "#FF61CC",
];

function getRegionColor(region) {
  if (regionColors[region]) return regionColors[region];
  if (region === "Other") return regionColors.Other;
  // cycle through palette if unknown region
  var keys = Object.keys(regionColors);
  var index = keys.length % ggplotPalette.length;
  var color = ggplotPalette[index];
  regionColors[region] = color;
  return color;
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
    h = s = 0;
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

function hslToHex(h, s, l) {
  h /= 360;
  const hue2rgb = (p, q, t) => {
    if (t < 0) t += 1;
    if (t > 1) t -= 1;
    if (t < 1 / 6) return p + (q - p) * 6 * t;
    if (t < 1 / 2) return q;
    if (t < 2 / 3) return p + (q - p) * (2 / 3 - t) * 6;
    return p;
  };
  let r, g, b;
  if (s === 0) {
    r = g = b = l;
  } else {
    const q = l < 0.5 ? l * (1 + s) : l + s - l * s;
    const p = 2 * l - q;
    r = hue2rgb(p, q, h + 1 / 3);
    g = hue2rgb(p, q, h);
    b = hue2rgb(p, q, h - 1 / 3);
  }
  const toHex = (x) =>
    Math.round(x * 255)
      .toString(16)
      .padStart(2, "0");
  return "#" + toHex(r) + toHex(g) + toHex(b);
}

function gradientShades(hex, steps = 4) {
  const hsl = hexToHsl(hex);
  const shades = [];
  for (let i = 0; i < steps; i++) {
    const factor = (i / Math.max(1, steps - 1)) * 0.5;
    const l = Math.min(1, Math.max(0, hsl.l + factor));
    shades.push(hslToHex(hsl.h, hsl.s, l));
  }
  return shades;
}

function randomHex() {
  return (
    "#" +
    Math.floor(Math.random() * 16777215)
      .toString(16)
      .padStart(6, "0")
  );
}
