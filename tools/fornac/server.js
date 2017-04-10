var jsdom = require('jsdom');
var fs = require("fs");
var d3 = fs.readFileSync("././d3.js", "utf-8");
var fornac = fs.readFileSync("./fornac.js", "utf-8");

var config = {};

config.file = './editor.html';
config.features = {
	FetchExternalResources : ['script'],
	ProcessExternalResources : ['script']
};

config.done = function (err, window) {
  	// node server.js > test.svg to export
  	// add styles to the styles element in the svg
  	var svgStyle = window.$("style").html();
  	var svg = window.$("svg").attr("viewBox", "0 0 300 300");
  	svg.attr("xmlns","http://www.w3.org/2000/svg");
  	svg.attr("version","1.1");
  	svg.attr("xmlns:xlink","http://www.w3.org/1999/xlink");
  	window.$("svg > style").html(''+svgStyle);
  	var svgFinal = 
  		`<?xml version="1.0" standalone="no"?>
<?xml-stylesheet type="text/css" href="style8.css"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
` +
  		window.$("#rna_ss").html();
    console.log(svgFinal);
}

jsdom.env(config);