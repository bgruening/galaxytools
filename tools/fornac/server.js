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
  	// node service.js > test.svg to export
  	// add styles to the styles element in the svg

    console.log(window.$("#rna_ss").html());
}

jsdom.env(config);