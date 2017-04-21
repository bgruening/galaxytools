// Usage node server.js [format] [input] [output]

if(process.argv.length < 5){
	console.log("Script expects 5 or more arguments.")
	console.log("Usage: node server.js [format] [input] [output]");
	process.exit();
}else{

	// Requirements
	var jsdom = require('jsdom');
	var fs = require("fs");
	const pnfs = require("pn/fs");
	const svg2png = require("svg2png");

	// Global variables
	var editorHtml;
	var tmpFileName;
	var format = process.argv[2];
	var input = process.argv[3];
	var output = process.argv[4];
	var tool_directory = process.argv[5];

	function readInput() {
		return fs.readFileSync(input, "utf-8");
	}

	/** A function that loads the input file and then modifies the editor.html
	*	to change the input parameters for Fornac. It loads the structures and
	*	the nucleotides for the visualization.
	*/
	function writeInput(content) {
		content = JSON.stringify(content);
		editorHtml = fs.readFileSync(tool_directory + "editor.html", "utf-8");
		var newValue = editorHtml.replace(/<input>/, content);
		var removedExtension = cleanFileName(input);
		tmpFileName = tool_directory + '/' + 'editor_' + removedExtension + '.html';
	  	fs.writeFileSync(tmpFileName, newValue, 'utf-8');
	}

	/** Removal of the prefixes of the full path and extension trimming.
	*/
	function cleanFileName(name){
		var parts = name.split("\/");
		var nameAndExtension = parts[parts.length - 1];
		return nameAndExtension.split("\.")[0];
	}

	/** A function that deletes the temporary file that is used for the 
	*	particular input sequence. This file is generated from the writeInput()
	*	function.
	*/
	function deleteTmpFile() {
		fs.unlink(tmpFileName);
	}

	/** Writes the output of the final file. Depending on the format it might be
	*	an svg, png or some other format file.
	*/
	function writeSvgOutput(content){
		fs.writeFileSync(output, content, 'utf-8');
	}

	/** Processes and generates the svg, also adds the needed tags for it to be
	*	displayable as a xml file in browser.
	*/
	function generateSvg(window){
		// Get the css styles for the svg
	  	var svgStyle = window.$("style").html();
	  	// Svg attributes that are needed for standalone function
	  	// as a separate file
	  	var svg = window.$("svg").attr("viewBox", "0 0 300 300");
	  	svg.attr("xmlns","http://www.w3.org/2000/svg");
	  	svg.attr("version","1.1");
	  	svg.attr("xmlns:xlink","http://www.w3.org/1999/xlink");
	  	svg.attr("width","300");
	  	svg.attr("height","300");
	  	// Locate the style inside the svg itself and fill it with
	  	// the previously loaded styles
	  	window.$("svg > style").html(''+svgStyle);
	  	var svgFinal = 
	  		`<?xml version="1.0" standalone="no"?>
			<?xml-stylesheet type="text/css" href="style8.css"?>
			<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
			` +
  			window.$("#rna_ss").html();
	    writeSvgOutput(svgFinal);
	}

	/** Processes and generates the svg, also adds the needed tags for it to be
	*	displayable as a xml file in browser.
	*/
	function convertSvgToPng(){
		pnfs.readFile(output + ".svg")
		    .then(svg2png)
		    .then(buffer => fs.writeFile(output + ".png", buffer))
		    .catch(e => console.error(e));
	}

	/** A function that runs the JSDOM server and executes it with the needed options.
	*	The server is set to execute the javascripts so the svg can be generated.
	*/
	function runJSDOM(content){
		writeInput(content);

		// The config object is built for the jsdom options
		var config = {};

		config.file = tmpFileName;
		config.features = {
			FetchExternalResources : ['script'],
			ProcessExternalResources : ['script']
		};

		config.done = function (err, window) {
		  	if(format == "svg"){
		  		generateSvg(window);
			}else if(format == "png"){
			 	generateSvg(window);
			 	convertSvgToPng();
				fs.unlink(output + ".svg");
			}
			deleteTmpFile();
		}

		// Run the jsdom server and generate the output
		jsdom.env(config);
	}

	/** The function that contains the main logic.
	*/
	function main(){
		var content = JSON.parse(readInput());

		if(content.constructor === Array){
			for(object in content){
				runJSDOM(content);
			}
		}else{
			runJSDOM(content);
		}
	}

	// Call to main
	main();
}