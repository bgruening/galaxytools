// Usage node server.js [format] [input] [output]

if(process.argv.length < 5){
	console.log("Script expects 5 or more arguments.")
	console.log("Usage: node server.js [format] [input] [output]");
	process.exit();
}else{

	// requirements
	var jsdom = require('jsdom');
	var fs = require("fs");

	// global variables
	var editorHtml;
	var format = process.argv[2];
	var input = process.argv[3];
	var output = process.argv[4];

	/** A function that loads the input file and then modifies the editor.html
	*	to change the input parameters for Fornac. It loads the structures and
	*	the nucleotides for the visualization.
	*/
	function writeInput() {
		var inputFile = fs.readFileSync(input.toString(), "utf-8");
		editorHtml = fs.readFileSync("./editor.html", "utf-8");
		var newValue = editorHtml.replace(/<input>/, inputFile);

	  	fs.writeFileSync('./editor.html', newValue, 'utf-8');
	}

	/** Reverts the changes made to editor.html for future use. If this function
	*	is not applied the file will have to be changed manually for the script
	*	to work.
	*/
	function revertInput() {
		fs.writeFileSync('./editor.html', editorHtml, 'utf-8');
	}

	/** Writes the output of the final file. Depending on the format it might be
	*	an svg, png or some other format file.
	*/
	function writeOutput(content){
		fs.writeFileSync(output.toString(), content, 'utf-8');
	}

	/** The function that contains the main logic
	*/
	function main(){
		writeInput();

		// The config object is built for the jsdom options
		var config = {};

		config.file = './editor.html';
		config.features = {
			FetchExternalResources : ['script'],
			ProcessExternalResources : ['script']
		};

		config.done = function (err, window) {
		  	if(format == "svg"){
		  		// Get the css styles for the svg
			  	var svgStyle = window.$("style").html();
			  	// Svg attributes that are needed for standalone function
			  	// as a separate file
			  	var svg = window.$("svg").attr("viewBox", "0 0 300 300");
			  	svg.attr("xmlns","http://www.w3.org/2000/svg");
			  	svg.attr("version","1.1");
			  	svg.attr("xmlns:xlink","http://www.w3.org/1999/xlink");
			  	// Locate the style inside the svg itself and fill it with
			  	// the previously loaded styles
			  	window.$("svg > style").html(''+svgStyle);

			  	var svgFinal = 
			  		`<?xml version="1.0" standalone="no"?>
					<?xml-stylesheet type="text/css" href="style8.css"?>
					<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
					` +
		  			window.$("#rna_ss").html();
			    writeOutput(svgFinal);
			}else if(format == "png"){
				//TODO: handle png format request
			}
			revertInput();
		}

		// Run the jsdom server and generate the output
		jsdom.env(config);
	}

	// call to main
	main();
}