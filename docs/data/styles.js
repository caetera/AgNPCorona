var styles = [ {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.6.0",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "pH experiment",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "font-family" : "SansSerif.plain",
      "font-weight" : "normal",
      "shape" : "ellipse",
      "height" : 50.0,
      "background-opacity" : 1.0,
      "border-width" : 0.0,
      "border-color" : "rgb(204,204,204)",
      "text-valign" : "center",
      "text-halign" : "center",
      "color" : "rgb(0,0,0)",
      "border-opacity" : 1.0,
      "background-color" : "rgb(153,153,153)",
      "text-opacity" : 1.0,
      "width" : 50.0,
      "font-size" : 7,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[Class = 'Composition']",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[Class = 'Other properties']",
    "css" : {
      "shape" : "octagon"
    }
  }, {
    "selector" : "node[Class = 'Hydrophobicity']",
    "css" : {
      "shape" : "vee"
    }
  }, {
    "selector" : "node[Class = 'Alpha and turn propensities']",
    "css" : {
      "shape" : "diamond"
    }
  }, {
    "selector" : "node[Class = 'Beta propensity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[Class = 'Physicochemical properties']",
    "css" : {
      "shape" : "triangle"
    }
  }, {
    "selector" : "node[Color_pH > 1.58243101]",
    "css" : {
      "background-color" : "rgb(178,24,43)"
    }
  }, {
    "selector" : "node[Color_pH = 1.58243101]",
    "css" : {
      "background-color" : "rgb(255,0,0)"
    }
  }, {
    "selector" : "node[Color_pH > 1][Color_pH < 1.58243101]",
    "css" : {
      "background-color" : "mapData(Color_pH,1,1.58243101,rgb(153,153,153),rgb(255,0,0))"
    }
  }, {
    "selector" : "node[Color_pH > 0.62535073][Color_pH < 1]",
    "css" : {
      "background-color" : "mapData(Color_pH,0.62535073,1,rgb(0,255,0),rgb(153,153,153))"
    }
  }, {
    "selector" : "node[Color_pH = 0.62535073]",
    "css" : {
      "background-color" : "rgb(0,255,0)"
    }
  }, {
    "selector" : "node[Color_pH < 0.62535073]",
    "css" : {
      "background-color" : "rgb(247,247,247)"
    }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "target-arrow-color" : "rgb(0,0,0)",
      "font-family" : "Dialog.plain",
      "font-weight" : "normal",
      "text-opacity" : 1.0,
      "line-color" : "rgb(132,132,132)",
      "line-style" : "solid",
      "opacity" : 1.0,
      "font-size" : 10,
      "source-arrow-shape" : "none",
      "width" : 2.0,
      "color" : "rgb(0,0,0)",
      "source-arrow-color" : "rgb(0,0,0)",
      "target-arrow-shape" : "none",
      "content" : ""
    }
  }, {
    "selector" : "edge[Correlation > 1]",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[Correlation = 1]",
    "css" : {
      "width" : 20.0
    }
  }, {
    "selector" : "edge[Correlation > 0.89999998][Correlation < 1]",
    "css" : {
      "width" : "mapData(Correlation,0.89999998,1,10.0,20.0)"
    }
  }, {
    "selector" : "edge[Correlation > 0.54559425][Correlation < 0.89999998]",
    "css" : {
      "width" : "mapData(Correlation,0.54559425,0.89999998,5.0,10.0)"
    }
  }, {
    "selector" : "edge[Correlation = 0.54559425]",
    "css" : {
      "width" : 5.0
    }
  }, {
    "selector" : "edge[Correlation < 0.54559425]",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.6.0",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Temperature experiment",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "font-family" : "SansSerif.plain",
      "font-weight" : "normal",
      "shape" : "ellipse",
      "height" : 50.0,
      "background-opacity" : 1.0,
      "border-width" : 0.0,
      "border-color" : "rgb(204,204,204)",
      "text-valign" : "center",
      "text-halign" : "center",
      "color" : "rgb(0,0,0)",
      "border-opacity" : 1.0,
      "background-color" : "rgb(153,153,153)",
      "text-opacity" : 1.0,
      "width" : 50.0,
      "font-size" : 7,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[Class = 'Undefined']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[Class = 'Composition']",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[Class = 'Other properties']",
    "css" : {
      "shape" : "octagon"
    }
  }, {
    "selector" : "node[Class = 'Hydrophobicity']",
    "css" : {
      "shape" : "vee"
    }
  }, {
    "selector" : "node[Class = 'Alpha and turn propensities']",
    "css" : {
      "shape" : "diamond"
    }
  }, {
    "selector" : "node[Class = 'Beta propensity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[Class = 'Physicochemical properties']",
    "css" : {
      "shape" : "triangle"
    }
  }, {
    "selector" : "node[Color_Temp > 1.47529729]",
    "css" : {
      "background-color" : "rgb(178,24,43)"
    }
  }, {
    "selector" : "node[Color_Temp = 1.47529729]",
    "css" : {
      "background-color" : "rgb(255,0,0)"
    }
  }, {
    "selector" : "node[Color_Temp > 1.00000001][Color_Temp < 1.47529729]",
    "css" : {
      "background-color" : "mapData(Color_Temp,1.00000001,1.47529729,rgb(153,153,153),rgb(255,0,0))"
    }
  }, {
    "selector" : "node[Color_Temp > 0.59464339][Color_Temp < 1.00000001]",
    "css" : {
      "background-color" : "mapData(Color_Temp,0.59464339,1.00000001,rgb(0,255,0),rgb(153,153,153))"
    }
  }, {
    "selector" : "node[Color_Temp = 0.59464339]",
    "css" : {
      "background-color" : "rgb(0,255,0)"
    }
  }, {
    "selector" : "node[Color_Temp < 0.59464339]",
    "css" : {
      "background-color" : "rgb(247,247,247)"
    }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "target-arrow-color" : "rgb(0,0,0)",
      "font-family" : "Dialog.plain",
      "font-weight" : "normal",
      "text-opacity" : 1.0,
      "line-color" : "rgb(132,132,132)",
      "line-style" : "solid",
      "opacity" : 1.0,
      "font-size" : 10,
      "source-arrow-shape" : "none",
      "width" : 2.0,
      "color" : "rgb(0,0,0)",
      "source-arrow-color" : "rgb(0,0,0)",
      "target-arrow-shape" : "none",
      "content" : ""
    }
  }, {
    "selector" : "edge[Correlation > 1]",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[Correlation = 1]",
    "css" : {
      "width" : 20.0
    }
  }, {
    "selector" : "edge[Correlation > 0.89999998][Correlation < 1]",
    "css" : {
      "width" : "mapData(Correlation,0.89999998,1,10.0,20.0)"
    }
  }, {
    "selector" : "edge[Correlation > 0.54559425][Correlation < 0.89999998]",
    "css" : {
      "width" : "mapData(Correlation,0.54559425,0.89999998,5.0,10.0)"
    }
  }, {
    "selector" : "edge[Correlation = 0.54559425]",
    "css" : {
      "width" : 5.0
    }
  }, {
    "selector" : "edge[Correlation < 0.54559425]",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
} ]
