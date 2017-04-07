/* fornai.js
* A container for display RNA secondary structure.
*
* Author: Peter Kerpedjiev <pkerp@tbi.univie.ac.at>
* Version: 0.2
* Date: 2015-03-15
*/

function FornaContainer(element, passedOptions) {
    var self = this;

    self.options = {
        "displayAllLinks": false,
        "labelInterval": 10,
        "applyForce": true,
        "initialSize": [200,200],
        "allowPanningAndZooming": true
    };

    if (arguments.length > 1) {
        for (var option in passedOptions) {
            if (self.options.hasOwnProperty(option))
                self.options[option] = passedOptions[option];
        }
    }

    self.options.svgW = self.options.initialSize[0];
    self.options.svgH = self.options.initialSize[1];

    var fill = d3.scale.category20();

    // mouse event vars
    var mousedownLink = null,
        mousedownNode = null,
        mouseupNode = null;

    var xScale = d3.scale.linear()
    .domain([0,self.options.svgW]).range([0,self.options.svgW]);
    var yScale = d3.scale.linear()
    .domain([0,self.options.svgH]).range([0, self.options.svgH]);

    var graph = self.graph = {
        "nodes":[],
        "links":[]
    };
    
    self.linkStrengths = {
        "pseudoknot": 0.5,
        "basepair": 0.5,
        "proteinChain": 0.00,
        "chainChain": 0.00,
        "intermolecule": 10.00,
        "other": 10.00
    };
    
    self.displayParameters = {
        "displayBackground": "true",
        "displayNumbering": "true",
        "displayNodeOutline": "true",
        "displayNodeLabel": "true",
        "displayLinks": "true",
        "displayPseudoknotLinks": "true",
        "displayProteinLinks": "true"
    };

    self.colorScheme = 'structure';
    self.customColors = {};
    self.animation = self.options.applyForce;
    // don't listen to events because a model window is open somewhere
    self.deaf = false;
    self.rnas = {};
    self.extraLinks = []; //store links between different RNAs


    self.createInitialLayout = function(structure, passedOptions) {
        // the default options
        var options = { 
                        'sequence': '',
                        'name': 'empty',
                        'positions': [],
                        'labelInterval': 10,
                        'avoidOthers': true
                      };

        if (arguments.length == 2) {
            for (var option in passedOptions) {
                if (options.hasOwnProperty(option))
                    options[option] = passedOptions[option];
            }
        }

        rg = new RNAGraph(options.sequence, structure, options.name);
        rnaJson = rg.recalculateElements()

        if (options.positions.length === 0) {
            // no provided positions means we need to calculate an initial layout
            options.positions = simpleXyCoordinates(rnaJson.pairtable);
        }

        rnaJson = rnaJson.elementsToJson()
        .addPositions("nucleotide", options.positions)
        .addLabels(1, options.labelInterval)
        .reinforceStems()
        .reinforceLoops()
        .connectFakeNodes()

        return rnaJson;
    }

    self.addRNA = function(structure, passedOptions) {
        var rnaJson = self.createInitialLayout(structure, passedOptions);

        if (arguments.length === 1)
            passedOptions = {}

        if ('avoidOthers' in passedOptions)
            self.addRNAJSON(rnaJson, passedOptions.avoidOthers);
        else
            self.addRNAJSON(rnaJson, true);

        return rnaJson;
    }

    self.addRNAJSON = function(rnaGraph, avoidOthers) {
        // Add an RNAGraph, which contains nodes and links as part of the
        // structure
        // Each RNA will have uid to identify it
        // when it is modified, it is replaced in the global list of RNAs
        //
        var maxX, minX;

        if (avoidOthers) {
            if (self.graph.nodes.length > 0)
                maxX = d3.max(self.graph.nodes.map(function(d) { return d.x; }));
            else
                maxX = 0;

            minX = d3.min(rnaGraph.nodes.map(function(d) { return d.x; })); 

            rnaGraph.nodes.forEach(function(node) {
                node.x += (maxX - minX);
                node.px += (maxX - minX);
            });
        }

        rnaGraph.nodes.forEach(function(node) {
            node.rna = rnaGraph;
        });

        self.rnas[rnaGraph.uid] = rnaGraph;
        self.recalculateGraph();

        self.update();
        self.centerView();
    };

    self.transitionRNA = function(previousRNAJson, newStructure, options) {
        //transition from an RNA which is already displayed to a new structure
        var newRNAJson = self.createInitialLayout(newStructure, options);
        console.log('newRNAJson:', newRNAJson);

        visNodes.selectAll('g.gnode').each(function(d) { console.log('d before', d); });
        var gnodes = visNodes.selectAll('g.gnode').data(newRNAJson);

        gnodes.each(function(d) { console.log('d after', d); });

        gnodes.transition().attr('transform', function(d) { 
            console.log('d after', d);
            return 'translate(' + [d.x, d.y] + ')'}).duration(1000);

    };

    self.recalculateGraph = function(rnaGraph) {
        // Condense all of the individual RNAs into one
        // collection of nodes and links
        self.graph.nodes = [];
        self.graph.links = [];
        for (var uid in self.rnas) {
            self.graph.nodes = self.graph.nodes.concat(self.rnas[uid].nodes);
            self.graph.links = self.graph.links.concat(self.rnas[uid].links);
        }

        // Create a lookup table so that we can access each node
        // based on its uid. This will be used to create the links
        // between different RNAs
        var uidsToNodes = {};

        for (var i = 0; i < self.graph.nodes.length; i++)
            uidsToNodes[self.graph.nodes[i].uid] = self.graph.nodes[i];

        self.graph.links.forEach(function(link) {
            link.source = uidsToNodes[link.source.uid];
            link.target = uidsToNodes[link.target.uid];
        });

        for (i = 0; i < self.extraLinks.length; i++) {
            // the actual node objects may have changed, so we hae to recreate
            // the extra links based on the uids

            if (!(self.extraLinks[i].target.uid in uidsToNodes)) {
                console.log("not there:", self.extraLinks[i]);
            }

            self.extraLinks[i].source = uidsToNodes[self.extraLinks[i].source.uid];
            self.extraLinks[i].target = uidsToNodes[self.extraLinks[i].target.uid];
            
            if (self.extraLinks[i].linkType == 'intermolecule') {
                //remove links to middle nodes
                fakeLinks = self.graph.links.filter(function(d) { 
                    return ((d.source == self.extraLinks[i].source || d.source == self.extraLinks[i].target ||
                            d.target == self.extraLinks[i].source || d.target == self.extraLinks[i].source) &&
                            d.linkType == 'fake');
                });

                for (var j = 0; j < fakeLinks.length; j++) {
                    var linkIndex = self.graph.links.indexOf(fakeLinks[j]); 
                    self.graph.links.splice(linkIndex, 1);
                }
            }

            graph.links.push(self.extraLinks[i]);
        }
    };

    self.addNodes = function addNodes(json) {
        // add a new set of nodes from a json file

        // Resolve the sources and targets of the links so that they
        // are not just indeces into an array
        json.links.forEach(function(entry) {
            if (typeof entry.source == "number") entry.source = json.nodes[entry.source];
            if (typeof entry.target == "number") entry.target = json.nodes[entry.target];
        });

        // Get the maximum x and y values of the current graph
        // so that we don't place a new structure on top of the
        // old one
        if (self.graph.nodes.length > 0) {
            maxX = d3.max(self.graph.nodes.map(function(d) {return d.x;}));
            maxY = d3.max(self.graph.nodes.map(function(d) {return d.y;}));
        } else {
            maxX = 0;
            maxY = 0;
        }

        json.nodes.forEach(function(entry) {
            if (!(entry.rna.uid in self.rnas)) {
                self.rnas[entry.rna.uid] = entry.rna;
            }

            entry.x += maxX;
            //entry.y += maxY;

            entry.px += maxX;
            //entry.py += maxY;
        });

        r = new RNAGraph('','');
        r.nodes = json.nodes;
        r.links = json.links;

        //self.addRNA(r);
        self.recalculateGraph();

        self.update();
        self.centerView();
    };

    self.addCustomColors = function addCustomColors(json) {
        // Add a json file containing the custom colors
        self.customColors = json;
    };

    self.clearNodes = function clearNodes() {
        self.graph.nodes = [];
        self.graph.links = [];

        self.rnas = {};
        self.extraLinks = [];

        self.update();
    };
    
    self.toJSON = function toJSON() {
       var data = {"rnas": self.rnas, "extraLinks": self.extraLinks};
            var dataString = JSON.stringify(data, function(key, value) {
            //remove circular references
            if (key == 'rna') {
                return;
            } else {
                return value;
            }
       }, "\t");
       return dataString;
    };

    self.fromJSON = function(jsonString) {
        try{
            var data = JSON.parse(jsonString);
            var rnas = data.rnas;
            var extraLinks = data.extraLinks;
        } catch(err) {
            throw err;
        }

        for (uid in rnas) {
            if (rnas[uid].type == 'rna') {
                r = new RNAGraph()

                r.seq = rnas[uid].seq;
                r.dotbracket = rnas[uid].dotbracket;
                r.circular = rnas[uid].circular;
                r.pairtable = rnas[uid].pairtable;
                r.uid = rnas[uid].uid;
                r.structName = rnas[uid].structName;
                r.nodes = rnas[uid].nodes;
                r.links = rnas[uid].links;
                r.rnaLength = rnas[uid].rnaLength;
                r.elements = rnas[uid].elements;
                r.nucsToNodes = rnas[uid].nucsToNodes;
                r.pseudoknotPairs = rnas[uid].pseudoknotPairs;
            } else {
                r = new ProteinGraph()
                r.size = rnas[uid].size;
                r.nodes = rnas[uid].nodes;
                r.uid = rnas[uid].uid;
            }

            self.addRNAJSON(r, false);
        }

        extraLinks.forEach(function(link) {
            self.extraLinks.push(link);
        });

        self.recalculateGraph();
        self.update();
    };

    function setSize() {
        var svgW = $(element).width();
        var svgH = $(element).height();

        self.options.svgW = svgW;
        self.options.svgH = svgH;

        //Set the output range of the scales
        xScale.range([0, svgW]).domain([0, svgW]);
        yScale.range([0, svgH]).domain([0, svgH]);

        //re-attach the scales to the zoom behaviour
        zoomer.x(xScale)
        .y(yScale);

        self.brusher.x(xScale)
        .y(yScale);

        //resize the background
        rect.attr("width", svgW)
        .attr("height", svgH);

        svg.attr("width", svgW)
        .attr("height", svgH);

        self.centerView();
    }

    function changeColors(moleculeColors, d, scale) {
        if (moleculeColors.hasOwnProperty(d.num)) {
            val = parseFloat(moleculeColors[d.num]);

            if (isNaN(val)) {
                // passed in color is not a scalar, so 
                // treat it as a color
                return moleculeColors[d.num];
            } else {
                // the user passed in a float, let's use a colormap
                // to convert it to a color
                return scale(val);
            }
        } else {
            return 'white';
        }
    }

    self.changeColorScheme = function(newColorScheme) {
        var proteinNodes = visNodes.selectAll('[node_type=protein]');

        proteinNodes.classed("protein", true)
                    .attr('r', function(d) { return d.radius; });

        var gnodes = visNodes.selectAll('g.gnode');
        var circles = visNodes.selectAll('g.gnode').selectAll('circle');
        var nodes = visNodes.selectAll('g.gnode').select('[node_type=nucleotide]');
        self.colorScheme = newColorScheme;


        if (newColorScheme == 'sequence') {
            scale = d3.scale.ordinal()
            .range(['#dbdb8d', '#98df8a', '#ff9896', '#aec7e8', '#aec7e8'])
            .domain(['A','C','G','U','T']);
            nodes.style('fill', function(d) { 
                return scale(d.name);
            });

        } else if (newColorScheme == "structure") {
            scale = d3.scale.category10()
            .domain(['s','m','i','e','t','h','x'])
            .range(['lightgreen', '#ff9896', '#dbdb8d', 'lightsalmon',
                   'lightcyan', 'lightblue', 'transparent']);

                   nodes.style('fill', function(d) { 
                       return scale(d.elemType);
                   });

        } else if (newColorScheme == 'positions') {
            nodes.style('fill', function(d) { 
                scale = d3.scale.linear()
                .range(["#98df8a", "#dbdb8d", "#ff9896"])
                .interpolate(d3.interpolateLab)
                .domain([1, 1 + (d.rna.rnaLength - 1) / 2, d.rna.rnaLength]);

                return scale(d.num);
            });
        } else if (newColorScheme == 'custom') {
            // scale to be used in case the user passes scalar
            // values rather than color names
            scale = d3.scale.linear()
            .interpolate(d3.interpolateLab)
            .domain(self.customColors.domain)
            .range(self.customColors.range);

            nodes.style('fill', function(d) {
                if (typeof self.customColors == 'undefined') {
                    return 'white';
                }
                
                if (self.customColors.colorValues.hasOwnProperty(d.structName) &&
                    self.customColors.colorValues[d.structName].hasOwnProperty(d.num)) {
                    // if a molecule name is specified, it supercedes the default colors
                    // (for which no molecule name has been specified)
                    moleculeColors = self.customColors.colorValues[d.structName];
                    return changeColors(moleculeColors, d, scale);
                } else if (self.customColors.colorValues.hasOwnProperty('')) {
                    moleculeColors = self.customColors.colorValues[''];
                    return changeColors(moleculeColors, d, scale);
                }

                return 'white';
            });
        }
    };

    function mousedown() {

    }

    function mousemove() {
        if (!mousedownNode) return;

        mpos = d3.mouse(vis.node());
        // update drag line
        dragLine
        .attr("x1", mousedownNode.x)
        .attr("y1", mousedownNode.y)
        .attr("x2", mpos[0])
        .attr("y2", mpos[1]);

    }

    function mouseup() {
        if (mousedownNode) {
            dragLine
            .attr("class", "dragLineHidden");
        }

        // clear mouse event vars
        resetMouseVars();
        //update()
    }
    //adapt size to window changes:
    window.addEventListener("resize", setSize, false);

    zoomer = d3.behavior.zoom()
        .scaleExtent([0.1,10])
        .x(xScale)
        .y(yScale)
        .on("zoomstart", zoomstart)
        .on("zoom", redraw);

    d3.select(element).select("svg").remove();

    var svg = d3.select(element)
    .attr("tabindex", 1)
    .on("keydown.brush", keydown)
    .on("keyup.brush", keyup)
    .each(function() { this.focus(); })
    .append("svg:svg")
    .attr("width", self.options.svgW)
    .attr("height", self.options.svgH)
    .attr("id", 'plotting-area');

    // set css for svg
    var style = svg.append('svg:style');
    $.get(location.origin + window.location.pathname.replace(/[^\/]+$/g,"") + "css/fornac.css", function(content){
        style.text(content.replace(/[\s\n]/g, ""));
    });
    
    self.options.svg = svg;

    var svgGraph = svg.append('svg:g')
    .on('mousemove', mousemove)
    .on('mousedown', mousedown)
    .on('mouseup', mouseup)

    if (self.options.allowPanningAndZooming)
        svgGraph.call(zoomer)

    var rect = svgGraph.append('svg:rect')
    .attr('width', self.options.svgW)
    .attr('height', self.options.svgH)
    .attr('fill', 'white')
    .attr('stroke', 'transparent')
    .attr('stroke-width', 1)
    //.attr("pointer-events", "all")
    .attr("id", "zrect");

    var brush = svgGraph.append('g')
    .datum(function() { return {selected: false, previouslySelected: false}; })
    .attr("class", "brush");
    var vis = svgGraph.append("svg:g");
    var visLinks = vis.append("svg:g");
    var visNodes = vis.append("svg:g");


    self.brusher = d3.svg.brush()
                .x(xScale)
                .y(yScale)
               .on("brushstart", function(d) {
                   var gnodes = visNodes.selectAll('g.gnode').selectAll('.outlineNode');
                   gnodes.each(function(d) { d.previouslySelected = ctrlKeydown && d.selected; });
               })
               .on("brush", function() {
                   var gnodes = visNodes.selectAll('g.gnode').selectAll('.outlineNode');
                   var extent = d3.event.target.extent();

                   gnodes.classed("selected", function(d) {
                       return d.selected = self.options.applyForce && d.previouslySelected ^
                       (extent[0][0] <= d.x && d.x < extent[1][0]
                        && extent[0][1] <= d.y && d.y < extent[1][1]);
                   });
               })
               .on("brushend", function() {
                   d3.event.target.clear();
                   d3.select(this).call(d3.event.target);
               })

      brush.call(self.brusher)
          .on("mousedown.brush", null)
          .on("touchstart.brush", null)                                                                      
          .on("touchmove.brush", null)                                                                       
          .on("touchend.brush", null);                                                                       
      brush.select('.background').style('cursor', 'auto')

    function zoomstart() {
        var node = visNodes.selectAll('g.gnode').selectAll('.outlineNode');
        node.each(function(d) {
                d.selected = false;
                d.previouslySelected = false;
                })
        node.classed("selected", false);
    }

    function redraw() {
        vis.attr("transform",
                 "translate(" + d3.event.translate + ")" + " scale(" + d3.event.scale + ")");
    }

    self.centerView = function() {
        // Center the view on the molecule(s) and scale it so that everything
        // fits in the window

        //no molecules, nothing to do
        if (self.graph.nodes.length === 0)
            return;

        // Get the bounding box
        minX = d3.min(self.graph.nodes.map(function(d) {return d.x;}));
        minY = d3.min(self.graph.nodes.map(function(d) {return d.y;}));

        maxX = d3.max(self.graph.nodes.map(function(d) {return d.x;}));
        maxY = d3.max(self.graph.nodes.map(function(d) {return d.y;}));


        // The width and the height of the molecule
        molWidth = maxX - minX;
        molHeight = maxY - minY;

        // how much larger the drawing area is than the width and the height
        widthRatio = self.options.svgW / (molWidth + 1);
        heightRatio = self.options.svgH / (molHeight + 1);

        // we need to fit it in both directions, so we scale according to
        // the direction in which we need to shrink the most
        minRatio = Math.min(widthRatio, heightRatio) * 0.8;

        // the new dimensions of the molecule
        newMolWidth = molWidth * minRatio;
        newMolHeight = molHeight * minRatio;

        // translate so that it's in the center of the window
        xTrans = -(minX) * minRatio + (self.options.svgW - newMolWidth) / 2;
        yTrans = -(minY) * minRatio + (self.options.svgH - newMolHeight) / 2;


        // do the actual moving
        vis.attr("transform",
                 "translate(" + [xTrans, yTrans] + ")" + " scale(" + minRatio + ")");

        // tell the zoomer what we did so that next we zoom, it uses the
        // transformation we entered here
        zoomer.translate([xTrans, yTrans ]);
        zoomer.scale(minRatio);

    };

    self.force = d3.layout.force()
    .charge(function(d) { if (d.nodeType == 'middle')  {
            return -30; 
    }
        else 
            return -30;})
    .chargeDistance(300)
    .friction(0.35)
    .linkDistance(function(d) { return 15 * d.value; })
    .linkStrength(function(d) { if (d.linkType in self.linkStrengths) {
                                  return self.linkStrengths[d.linkType];
                                } else {
                                  return self.linkStrengths.other; }
    })
    .gravity(0.000)
    .nodes(self.graph.nodes)
    .links(self.graph.links)
    .chargeDistance(110)
    .size([self.options.svgW, self.options.svgH]);

    // line displayed when dragging new nodes
    var dragLine = vis.append("line")
    .attr("class", "dragLine")
    .attr("x1", 0)
    .attr("y1", 0)
    .attr("x2", 0)
    .attr("y2", 0);

    function resetMouseVars() {
        mousedownNode = null;
        mouseupNode = null;
        mousedownLink = null;
    }

    var shiftKeydown = false;
    var ctrlKeydown = false;

    function selectedNodes(mouseDownNode) {
        var gnodes = visNodes.selectAll('g.gnode');

        if (ctrlKeydown) {
            return gnodes.filter(function(d) { return d.selected; });

            //return d3.selectAll('[structName=' + mouseDownNode.structName + ']');
        } else {
            return gnodes.filter(function(d) { return d.selected ; });
            //return d3.select(this);
        }
    }

    function dragstarted(d) {
        d3.event.sourceEvent.stopPropagation();

      if (!d.selected && !ctrlKeydown) {
          // if this node isn't selected, then we have to unselect every other node
            var node = visNodes.selectAll('g.gnode').selectAll('.outline_node');
            node.classed("selected", function(p) { return p.selected =  self.options.applyForce && (p.previouslySelected = false); })
          }

        d3.select(this).select('.outline_node').classed("selected", function(p) { d.previouslySelected = d.selected; return d.selected = self.options.applyForce && true; });

        var toDrag = selectedNodes(d);
        toDrag.each(function(d1) {
            d1.fixed |= 2;
        });

        //d3.event.sourceEvent.stopPropagation();
        //d3.select(self).classed("dragging", true);
        //
    }

    function dragged(d) {

        var toDrag = selectedNodes(d);

        toDrag.each(function(d1) {
            d1.x += d3.event.dx;
            d1.y += d3.event.dy;

            d1.px += d3.event.dx;
            d1.py += d3.event.dy;
        });

        self.resumeForce();
        d3.event.sourceEvent.preventDefault();
    }

    self.resumeForce = function() {
        if (self.animation)
            self.force.resume();
    }

    function dragended(d) {
        var toDrag = selectedNodes(d);

        toDrag.each(function(d1) {
            d1.fixed &= ~6;
        });
    }

    function collide(node) {
        var r = node.radius + 16,
        nx1 = node.x - r,
        nx2 = node.x + r,
        ny1 = node.y - r,
        ny2 = node.y + r;
        return function(quad, x1, y1, x2, y2) {
            if (quad.point && (quad.point !== node)) {
                var x = node.x - quad.point.x,
                y = node.y - quad.point.y,
                l = Math.sqrt(x * x + y * y),
                r = node.radius + quad.point.radius;
                if (l < r) {
                    l = (l - r) / l * .1;
                    node.x -= x *= l;
                    node.y -= y *= l;
                    quad.point.x += x;
                    quad.point.y += y;
                }
            }
            return x1 > nx2 || x2 < nx1 || y1 > ny2 || y2 < ny1;
        };
    }


    var drag = d3.behavior.drag()
    //.origin(function(d) { return d; })
    .on("dragstart", dragstarted)
    .on("drag", dragged)
    .on("dragend", dragended);

    function keydown() {
        if (self.deaf)
            // lalalalal, not listening
            return;

        if (shiftKeydown) return;

        keyIsDown = true;
        switch (d3.event.keyCode) {
            case 16:
                shiftKeydown = true;
                break;
            case 17:
                ctrlKeydown = true;
                break;
            case 67: //c
                self.centerView();
                break;
        }

        if (shiftKeydown || ctrlKeydown) {
            svgGraph.call(zoomer)
            .on("mousedown.zoom", null)
            .on("touchstart.zoom", null)
            .on("touchmove.zoom", null)
            .on("touchend.zoom", null);

            //svgGraph.on('zoom', null);
            vis.selectAll('g.gnode')
            .on('mousedown.drag', null);
        }

        if (ctrlKeydown) {
          brush.select('.background').style('cursor', 'crosshair')
          brush.call(self.brusher);
        }
    }

    function keyup() {
        shiftKeydown = false;
        ctrlKeydown = false;

        brush.call(self.brusher)
        .on("mousedown.brush", null)
        .on("touchstart.brush", null)                                                                      
        .on("touchmove.brush", null)                                                                       
        .on("touchend.brush", null);                                                                       

        brush.select('.background').style('cursor', 'auto')
        svgGraph.call(zoomer);

        vis.selectAll('g.gnode')
        .call(drag);
    }

    d3.select(element)
    .on('keydown', keydown)
    .on('keyup', keyup)
    .on('contextmenu', function() {
            d3.event.preventDefault(); 
    });

    linkKey = function(d) {
        return d.uid;
    };

    nodeKey = function(d) {
        key = d.uid;
        return key;
    };

    updateRnaGraph = function(r) {
        var nucleotidePositions = r.getPositions('nucleotide');
        var labelPositions = r.getPositions('label');

        var uids = r.getUids();

        r.recalculateElements()
        .elementsToJson()
        .addPseudoknots()
        .addPositions('nucleotide', nucleotidePositions)
        .addUids(uids)
        .addLabels(1, self.options.labelInterval)
        .addPositions('label', labelPositions)
        .reinforceStems()
        .reinforceLoops()
        .updateLinkUids()
    };

    removeLink = function(d) {
        // remove a link between two nodes
        index = self.graph.links.indexOf(d);

        if (index > -1) {
            //remove a link
            //graph.links.splice(index, 1);

            // there should be two cases
            // 1. The link is within a single molecule

            if (d.source.rna == d.target.rna) {
                var r = d.source.rna;

                r.addPseudoknots();
                r.pairtable[d.source.num] = 0;
                r.pairtable[d.target.num] = 0;

                updateRnaGraph(r);

            } else {
                // 2. The link is between two different molecules
                extraLinkIndex = self.extraLinks.indexOf(d);

                self.extraLinks.splice(extraLinkIndex, 1);
            }

            self.recalculateGraph();
        }

        self.update();
    };

    linkClick = function(d) {
        if (!shiftKeydown) {
            return;
        }

        var invalidLinks = {'backbone': true,
                             'fake': true,
                             'fakeFake': true,
                             'labelLink': true}

        if (d.linkType in invalidLinks ) 
            return;

        removeLink(d);
    };


    self.addLink =  function(newLink) {
        // this means we have a new json, which means we have
        // to recalculate the structure and change the colors
        // appropriately
        //
        if (newLink.source.rna == newLink.target.rna) {
            r = newLink.source.rna;

            r.pairtable[newLink.source.num] = newLink.target.num;
            r.pairtable[newLink.target.num] = newLink.source.num;

            updateRnaGraph(r);

        } else {
            //Add an extra link
            newLink.linkType = 'intermolecule';
            self.extraLinks.push(newLink);
        }
        self.recalculateGraph();
        self.update();
    };

    nodeMouseclick = function(d) {
        if (d3.event.defaultPrevented) return;

        if (!ctrlKeydown) {
            //if the shift key isn't down, unselect everything
            var node = visNodes.selectAll('g.gnode').selectAll('.outline_node');
            node.classed("selected", function(p) { return p.selected =  self.options.applyForce && (p.previouslySelected = false); })
        }

        // always select this node
        d3.select(this).select('circle').classed("selected", d.selected = self.options.applyForce && !d.previouslySelected);
    }

    nodeMouseup = function(d) {
        if (mousedownNode) {
            mouseupNode = d;

            if (mouseupNode == mousedownNode) { resetMouseVars(); return; }
            var newLink = {source: mousedownNode, target: mouseupNode, linkType: 'basepair', value: 1, uid:generateUUID()};

            for (i = 0; i < self.graph.links.length; i++) {
                if ((self.graph.links[i].source == mousedownNode)  || 
                    (self.graph.links[i].target == mousedownNode) ||
                        (self.graph.links[i].source == mouseupNode) ||
                            (self.graph.links[i].target == mouseupNode)) {

                    if (self.graph.links[i].linkType == 'basepair' || self.graph.links[i].linkType == 'pseudoknot') {
                        return;
                    }
                }

                if (((self.graph.links[i].source == mouseupNode)  && 
                     (self.graph.links[i].target == mousedownNode)) ||
                         ((self.graph.links[i].source == mousedownNode)  && 
                          (self.graph.links[i].target == mouseupNode))) {
                    if (self.graph.links[i].linkType == 'backbone') {
                        return;
                    }
                }
            }

            if (mouseupNode.nodeType == 'middle' || mousedownNode.nodeType == 'middle' || mouseupNode.nodeType == 'label' || mousedownNode.nodeType == 'label')
                return;

            self.addLink(newLink);

        }
    };

    nodeMousedown = function(d) {
      if (!d.selected && !ctrlKeydown) {
          // if this node isn't selected, then we have to unselect every other node
            var node = visNodes.selectAll('g.gnode').selectAll('.outlineNode');
            node.classed("selected", function(p) { return p.selected =  p.previouslySelected = false; })
          }


          d3.select(this).classed("selected", function(p) { d.previouslySelected = d.selected; return d.selected = self.options.applyForce && true; });

        if (!shiftKeydown) {
            return;
        }

        mousedownNode = d;

        dragLine
        .attr("class", "drag_line")
        .attr("x1", mousedownNode.x)
        .attr("y1", mousedownNode.y)
        .attr("x2", mousedownNode.x)
        .attr("y2", mousedownNode.y);

        //gnodes.attr('pointer-events',  'none');

    };

    self.startAnimation = function() {
      self.animation = true;
      vis.selectAll('g.gnode')
        .call(drag);
      self.force.start();
    };
    
    self.stopAnimation = function() {
      self.animation = false;
      vis.selectAll('g.gnode')
           .on('mousedown.drag', null);
      self.force.stop();
    };
    
    self.setFriction = function(value) {
      self.force.friction(value);
      self.resumeForce();
    };

    self.setCharge = function(value) {
      self.force.charge(value);
      self.resumeForce();
    };
    
    self.setGravity = function(value) {
      self.force.gravity(value);
      self.resumeForce();
    };
    
    self.setPseudoknotStrength = function(value) {
      self.linkStrengths.pseudoknot = value;
      self.update();
    };
    
    self.displayBackground = function(value) {
      self.displayParameters.displayBackground = value;
      self.updateStyle();
    };
    
    self.displayNumbering = function(value) {
      self.displayParameters.displayNumbering = value;
      self.updateStyle();
    };

    self.displayNodeOutline = function(value) {
      self.displayParameters.displayNodeOutline = value;
      self.updateStyle();
    };
    
    self.displayNodeLabel = function(value) {
      self.displayParameters.displayNodeLabel = value;
      self.updateStyle();
    };
    
    self.displayLinks = function(value) {
      self.displayParameters.displayLinks = value;
      self.updateStyle();
    };

    self.displayPseudoknotLinks = function(value) {
      self.displayParameters.displayPseudoknotLinks = value;
      self.updateStyle();
    };

    self.displayProteinLinks = function(value) {
      self.displayParameters.displayProteinLinks = value;
      self.updateStyle();
    };
    
    self.updateStyle = function() {
        // Background
        rect.classed("transparent", !self.displayParameters.displayBackground);
        // Numbering
        visNodes.selectAll('[node_type=label]').classed("transparent", !self.displayParameters.displayNumbering);
        visNodes.selectAll('[label_type=label]').classed("transparent", !self.displayParameters.displayNumbering);
        visLinks.selectAll('[link_type=label_link]').classed("transparent", !self.displayParameters.displayNumbering);
        // Node Outline
        svg.selectAll('circle').classed("hidden_outline", !self.displayParameters.displayNodeOutline);
        // Node Labels
        visNodes.selectAll('[label_type=nucleotide]').classed("transparent", !self.displayParameters.displayNodeLabel);
        // Links
        svg.selectAll("[link_type=real],[link_type=basepair],[link_type=backbone],[link_type=pseudoknot],[link_type=protein_chain],[link_type=chain_chain]").classed("transparent", !self.displayParameters.displayLinks);
        // Pseudoknot Links
        svg.selectAll("[link_type=pseudoknot]").classed("transparent", !self.displayParameters.displayPseudoknotLinks);
        // Protein Links
        svg.selectAll("[link_type=protein_chain]").classed("transparent", !self.displayParameters.displayProteinLinks);
        // Fake Links
        visLinks.selectAll("[link_type=fake]").classed("transparent", !self.options.displayAllLinks);
    };

    function nudge(dx, dy) {
        node.filter(function(d) { return d.selected; })
        .attr("cx", function(d) { return d.x += dx; })
        .attr("cy", function(d) { return d.y += dy; })

        link.filter(function(d) { return d.source.selected; })
        .attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; });

        link.filter(function(d) { return d.target.selected; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

        d3.event.preventDefault();
    }

    self.update = function () {
        self.force.nodes(self.graph.nodes)
        .links(self.graph.links);
        
        if (self.animation) {
          self.force.start();
        }

        var allLinks = visLinks.selectAll("line.link")
        .data(self.graph.links, linkKey);

        var linksEnter = allLinks.enter();

        linkLines = linksEnter.append("svg:line");

        linkLines.append("svg:title")
        .text(linkKey);

        linkLines
        .classed("link", true)
        .attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; })
        .attr("link_type", function(d) { return d.linkType; } )
        .attr("class", function(d) { return d3.select(this).attr('class') + " " + d.linkType; })
        .attr('pointer-events', function(d) { if (d.linkType == 'fake') return 'none'; else return 'all';});

        allLinks.attr('class', '')
        .classed('link', true)
        .attr("link_type", function(d) { return d.linkType; } )
        .attr("class", function(d) { return d3.select(this).attr('class') + " " + d.linkType; })

            allLinks.exit().remove();

            if (self.displayFakeLinks)
                xlink = allLinks;
            else
                xlink = visLinks.selectAll("[link_type=real],[link_type=pseudoknot],[link_type=protein_chain],[link_type=chain_chain],[link_type=label_link],[link_type=backbone],[link_type=basepair],[link_type=fake],[link_type=intermolecule]");

            domain = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
            var colors = d3.scale.category10().domain(domain);

            var gnodes = visNodes.selectAll('g.gnode')
            .data(self.graph.nodes, nodeKey);
            //.attr('pointer-events', 'all');

            gnodesEnter = gnodes.enter()
            .append('g')
            .classed('noselect', true)
            .classed('gnode', true)
             .attr('struct_name', function(d) { return d.structName; })
              .attr("transform", function(d) { 
                  if (typeof d.x != 'undefined' && typeof d.y != 'undefined')
                    return 'translate(' + [d.x, d.y] + ')'; 
                else
                    return ''
                })
             .each( function(d) { d.selected = d.previouslySelected = false; })

            gnodesEnter
            .call(drag)
            .on('mousedown', nodeMousedown)
            .on('mousedrag', function(d) {})
            .on('mouseup', nodeMouseup)
            .on('click', nodeMouseclick)
            .transition()
            .duration(750)
            .ease("elastic")
            .attr("r", 6.5);

            nodeTooltip = function(d) {
                nodeTooltips = {};

                nodeTooltips.nucleotide = d.num;
                nodeTooltips.label = '';
                nodeTooltips.pseudo = '';
                nodeTooltips.middle = '';
                nodeTooltips.protein = d.structName;

                return nodeTooltips[d.nodeType];
            };

            xlink.on('click', linkClick);

            var circleUpdate = gnodes.select('circle');

            // create nodes behind the circles which will serve to highlight them
            var nucleotideNodes = gnodesEnter.filter(function(d) { 
                return d.nodeType == 'nucleotide' || d.nodeType == 'label' || d.nodeType == 'protein';
            })
            nucleotideNodes.append("svg:circle")
            .attr('class', "outline_node")
            .attr("r", function(d) { return d.radius+1; })

            var node = gnodesEnter.append("svg:circle")
            .attr("class", "node")
            .classed("label", function(d) { return d.nodeType == 'label'; })
            .attr("r", function(d) { 
                if (d.nodeType == 'middle') return 0; 
                else {
                    return d.radius; 
                }
                })
            .attr("node_type", function(d) { return d.nodeType; })
            
            var labels = gnodesEnter.append("text")
            .text(function(d) { return d.name; })
            .attr('text-anchor', 'middle')
            .attr('font-size', 8.0)
            .attr('font-weight', 'bold')
            .attr('y', 2.5)
            .attr('class', 'node-label')
            .attr("label_type", function(d) { return d.nodeType; })
            .append("svg:title")
            .text(function(d) { 
                if (d.nodeType == 'nucleotide') {
                    return d.structName + ":" + d.num;
                } else {
                    return '';
                }
            });

            node.append("svg:title")
            .text(function(d) { 
                if (d.nodeType == 'nucleotide') {
                    return d.structName + ":" + d.num;
                } else {
                    return '';
                }
            });

            gnodes.exit().remove();

            //fake_nodes = self.graph.nodes.filter(function(d) { return d.nodeType == 'middle'; });
            //fake_nodes = self.graph.nodes.filter(function(d) { return true; });
            realNodes = self.graph.nodes.filter(function(d) { return d.nodeType == 'nucleotide' || d.nodeType == 'label';});

            self.force.on("tick", function() {

                var q = d3.geom.quadtree(realNodes),
                i = 0,
                n = realNodes.length;

                while (++i < n) q.visit(collide(realNodes[i]));

                xlink.attr("x1", function(d) { return d.source.x; })
                .attr("y1", function(d) {  return d.source.y; })
                .attr("x2", function(d) { return d.target.x; })
                .attr("y2", function(d) { return d.target.y; });

                // Translate the groups
                gnodes.attr("transform", function(d) { 
                    return 'translate(' + [d.x, d.y] + ')'; 
                });
            });
            
        self.changeColorScheme(self.colorScheme);

        if (self.animation) {
          self.force.start();
        }
        
        self.updateStyle();
    };
    
    setSize();
}
var numberSort = function(a,b) { return a - b; };

function generateUUID(){                                                                                        
    /* Stack Overflow:                                                                                          
     * http://stackoverflow.com/a/8809472/899470                                                                
     */
    var d = new Date().getTime();
    var uuid = 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
        var r = (d + Math.random()*16)%16 | 0;
        d = Math.floor(d/16);
        return (c=='x' ? r : (r&0x3|0x8)).toString(16);                                                         
    });                                                                                                         

    return uuid;
}

function isNormalInteger(str) {
    //http://stackoverflow.com/a/10834843/899470
    return /^\+?(0|[1-9]\d*)$/.test(str);
}

if(typeof(String.prototype.trim) === "undefined")
    {
        String.prototype.trim = function() 
        {
            return String(this).replace(/^\s+|\s+$/g, '');
        };
    }

function ColorScheme(colorsText) {
    var self = this;
    self.colorsText = colorsText;

    self.parseRange = function(rangeText) {
        //parse a number range such as 1-10 or 3,7,9 or just 7
        var parts = rangeText.split(',')
        var nums = [];

        for (var i = 0; i < parts.length; i++) {
            //could be 1 or 10-11  or something like that
            var parts1 = parts[i].split('-');

            if (parts1.length == 1)
                nums.push(parseInt(parts1[0]));
            else if (parts1.length == 2) {
                var from = parseInt(parts1[0]);
                var to = parseInt(parts1[1]);

                // add each number in this range
                for (var j = from; j <= to; j++) 
                    nums.push(j)
            } else {
                console.log('Malformed range (too many dashes):', rangeText);
            }
        }

        return nums;
    }

    self.parseColorText = function(colorText) {
        /* Parse the text of an RNA color string. Instructions and description
         * of the format are given below.
         *
         * The return is a json double dictionary indexed first by the 
         * molecule name, then by the nucleotide. This is then applied
         * by force.js to the RNAs it is displaying. When no molecule
         * name is specified, the color is applied to all molecules*/
        var lines = colorText.split('\n');
        var currMolecule = '';
        var counter = 1;
        var colorsJson = {colorValues: {'':{}}, range:['white', 'steelblue']};
        var domainValues = [];


        for (var i = 0; i < lines.length; i++) {

            if (lines[i][0] == '>') {
                // new molecule
                currMolecule = lines[i].trim().slice(1);
                counter = 1;

                colorsJson.colorValues[currMolecule] = {};
                continue;
            }

            words = lines[i].trim().split(/[\s]+/);

            for (var j = 0; j < words.length; j++) {
                if (isNaN(words[j])) {
                    if (words[j].search("range") === 0) {
                        //there's a color scale in this entry
                        parts = words[j].split('=');
                        partsRight = parts[1].split(':')
                        colorsJson.range = [partsRight[0], partsRight[1]];
                        continue;
                    }

                    if (words[j].search("domain") == 0) {
                        //there's a color scale in this entry
                        parts = words[j].split('=');
                        partsRight = parts[1].split(':')
                        colorsJson.domain = [partsRight[0], partsRight[1]];
                        continue;
                    }

                    // it's not a number, should be a combination 
                    // of a number (nucleotide #) and a color
                    parts = words[j].split(':');
                    nums = self.parseRange(parts[0]);
                    color = parts[1]

                    for (var k = 0; k < nums.length; k++) {
                        if (isNaN(color)) {
                            colorsJson.colorValues[currMolecule][nums[k]] = color;
                        } else {
                            colorsJson.colorValues[currMolecule][nums[k]] = +color;
                            domainValues.push(Number(color));
                        }
                    }
                } else {
                    //it's a number, so we add it to the list of values
                    //seen for this molecule
                    colorsJson.colorValues[currMolecule][counter] = Number(words[j]);
                    counter += 1;

                    domainValues.push(Number(words[j]));
                }
            }
        }

        if (!('domain' in colorsJson))
            colorsJson.domain = [Math.min.apply(null, domainValues), Math.max.apply(null, domainValues)];

        self.colorsJson = colorsJson;

        return self;
    };

    self.normalizeColors = function() {
        /* 
         * Normalize the passed in values so that they range from
         * 0 to 1
         */
        var value;

        for (var moleculeName in self.colorsJson) {
            var minNum = Number.MAX_VALUE;
            var maxNum = Number.MIN_VALUE;

            // iterate once to find the min and max values;
            for (var resnum in self.colorsJson.colorValues[moleculeName]) {
                value = self.colorsJson.colorValues[moleculeName][resnum];
                if (typeof value == 'number') {
                    if (value < minNum)
                        minNum = value;
                    if (value > maxNum)
                        maxNum = value;
                }
            }

            // iterate again to normalize
            for (resnum in self.colorsJson.colorValues[moleculeName]) {
                value = self.colorsJson.colorValues[moleculeName][resnum];
                if (typeof value == 'number') {
                    self.colorsJson.colorValues[moleculeName][resnum] = (value - minNum ) / (maxNum - minNum);
                }
            }
        }

        return self;
    };

    self.parseColorText(self.colorsText);
}

function ProteinGraph(structName, size, uid) {
    var self = this;

    self.type = 'protein';
    self.size = size;
    self.nodes = [{'name': 'P',
                   'num': 1,
                   'radius': 3 *  Math.sqrt(size),
                   'rna': self,
                   'nodeType': 'protein',
                   'structName': structName,
                   'elemType': 'p',
                   'size': size,
                   'uid': generateUUID()}];

    self.links = [];
    self.uid = generateUUID();

    self.addUids = function(uids) {
        for (var i = 0; i < uids.length; i++)
            self.nodes[i].uid = uids[i];

        return self;
    };

    self.getUids = function() {
        /* Get the positions of each node so that they
         * can be passed to elementsToJson later
         */
        uids = [];
        for (var i = 0; i < self.dotbracket.length; i++)
            uids.push(self.nodes[i].uid);

        return uids;
    };

}

function RNAGraph(seq, dotbracket, structName) {
    var self = this;

    self.type = 'rna';
    self.circularizeExternal = false;

    if (arguments.length == 0) {
        self.seq = '';
        self.dotbracket = '';
        self.structName = '';
    } else {
        self.seq = seq;
        self.dotbracket = dotbracket;  //i.e. ..((..))..
        self.structName = structName;
    }

    self.circular = false;

    if (self.dotbracket.length > 0 && self.dotbracket[self.dotbracket.length-1] == '*') {
        //circular RNA
        self.dotbracket = self.dotbracket.slice(0, self.dotbracket.length-1);
        self.circular = true;
    }

    self.uid = generateUUID();
    self.rnaLength = self.dotbracket.length;

    self.elements = [];            //store the elements and the 
                                   //nucleotides they contain
    self.pseudoknotPairs = [];
    self.nucsToNodes = {};

    self.addUids = function(uids) {
        for (var i = 0; i < uids.length; i++)
            self.nodes[i].uid = uids[i];

        return self;
    };

    self.computePairtable = function() {
        self.pairtable = rnaUtilities.dotbracketToPairtable(self.dotbracket);
    };

    self.computePairtable();

    self.addPositions = function(nodeType, positions) {
        labelNodes = self.nodes.filter(function(d) { return d.nodeType == nodeType; });

        for  (var i = 0; i < labelNodes.length; i++) {
            labelNodes[i].x = positions[i][0];
            labelNodes[i].px = positions[i][0];
            labelNodes[i].y = positions[i][1];
            labelNodes[i].py = positions[i][1];
        }

        return self;
    };

    self.getPositions = function(nodeType) {
        positions = [];
        nucleotideNodes = self.nodes.filter(function(d) { return d.nodeType == nodeType; });

        for (var i = 0; i < nucleotideNodes.length; i++)
            positions.push([nucleotideNodes[i].x, nucleotideNodes[i].y]);

        return positions;
    };

    self.getUids = function() {
        /* Get the positions of each node so that they
         * can be passed to elementsToJson later
         */
        uids = [];
        for (var i = 0; i < self.dotbracket.length; i++)
            uids.push(self.nodes[i].uid);

        return uids;
    };

    self.reinforceStems = function() {
        pt = self.pairtable;
        relevantElements = elements.filter( function(d) {
            return d[0] == 's' && d[2].length >= 4;
        });

        for (var i = 0; i < relevantElements.length; i++) {
            allNucs = relevantElements[i][2];
            nucs = allNucs.slice(0, allNucs.length / 2);

            for (var j = 0; j < nucs.length-1; j++) {
                self.addFakeNode([nucs[j], nucs[j+1], pt[nucs[j+1]], pt[nucs[j]]]);
            }
        }

        return self;    
    };

    self.reinforceLoops = function() {
        /* 
         * Add a set of fake nodes to enforce the structure
         */
        var filterNucs = function(d) { 
            return d !== 0 && d <= self.dotbracket.length;
        };

        for (i=0; i < self.elements.length; i++) {
            if (self.elements[i][0] == 's' || (!self.circularizeExternal && self.elements[i][0] == 'e'))
                continue;

            var nucs = self.elements[i][2].filter(filterNucs);

            if (self.elements[i][0] == 'e') {
                var newNode1 = {'name': '',
                    'num': -1,
                    //'radius': 18 * radius -6,
                    'radius': 0,
                    'rna': self,
                    'nodeType': 'middle',
                    'elemType': 'f',
                    'nucs': [],
                    'x': self.nodes[self.rnaLength-1].x,
                    'y': self.nodes[self.rnaLength-1].y,
                    'px': self.nodes[self.rnaLength-1].px,
                    'py': self.nodes[self.rnaLength-1].py,
                    'uid': generateUUID() };
                var newNode2 = {'name': '',
                    'num': -1,
                    //'radius': 18 * radius -6,
                    'radius': 0,
                    'rna': self,
                    'nodeType': 'middle',
                    'elemType': 'f',
                    'nucs': [],
                    'x': self.nodes[0].x,
                    'y': self.nodes[0].y,
                    'px': self.nodes[0].px,
                    'py': self.nodes[0].py,
                    'uid': generateUUID() };

                    nucs.push(self.nodes.length+1);
                    nucs.push(self.nodes.length+2);
                    self.nodes.push(newNode1);
                    self.nodes.push(newNode2);
            }
            

            self.addFakeNode(nucs);
        }

        return self;
    };

    self.updateLinkUids = function() {
        for (var i = 0; i < self.links.length; i++) {
            self.links[i].uid = self.links[i].source.uid + self.links[i].target.uid;
        }

        return self;
    }

    self.addFakeNode = function(nucs) {
        var linkLength = 18; //make sure this is consistent with the value in force.js
        var nodeWidth = 6;
        var angle = (3.1415 * 2) / (2 * nucs.length);
        var radius =  linkLength / (2 * Math.tan(angle));

        var fakeNodeUid = ""

        for (var i = 0; i < nucs.length; i++)
            fakeNodeUid += self.nodes[nucs[i]-1].uid;

        newNode = {'name': '',
                         'num': -1,
                         //'radius': 18 * radius -6,
                         'radius': radius,
                         'rna': self,
                         'nodeType': 'middle',
                         'elemType': 'f',
                         'nucs': nucs,
                         'uid': fakeNodeUid };
        self.nodes.push(newNode);

        newX = 0;
        newY = 0;
        coordsCounted = 0;

        angle = (nucs.length - 2) * 3.14159 / (2 * nucs.length);
        radius = 0.5 / Math.cos(angle);

        for (j = 0; j < nucs.length; j++) {
            if (nucs[j] === 0 || nucs[j] > self.dotbracket.length)
                continue;
            

            //link to the center node
            self.links.push({'source': self.nodes[nucs[j] - 1],
                             'target': self.nodes[self.nodes.length-1],
                             'linkType': 'fake',
                             'value': radius,
                             'uid': generateUUID() });

            if (nucs.length > 4) {
                //link across the loop
                self.links.push({'source': self.nodes[nucs[j] - 1],
                                 'target': self.nodes[nucs[(j + Math.floor(nucs.length / 2)) % nucs.length] - 1],
                                 'linkType': 'fake',
                                 'value': radius * 2,
                                 'uid': generateUUID() });
            }

            ia = ((nucs.length - 2) * 3.14159) / nucs.length;
            c = 2 * Math.cos(3.14159 / 2 - ia / 2);
            //link to over-neighbor
            self.links.push({'source': self.nodes[nucs[j] - 1],
                             'target': self.nodes[nucs[(j + 2) % nucs.length] - 1],
                             'linkType': 'fake',
                             'value': c});

            // calculate the mean of the coordinats in this loop
            // and place the fake node there
            fromNode = self.nodes[nucs[j]-1];
            if ('x' in fromNode) {
                newX += fromNode.x;
                newY += fromNode.y;

                coordsCounted += 1;
            }
        }

        if (coordsCounted > 0) {
            // the nucleotides had set positions so we can calculate the position
            // of the fake node
            newNode.x = newX / coordsCounted;
            newNode.y = newY / coordsCounted;
            newNode.px = newNode.x;
            newNode.py = newNode.y;
        }

        return self;

    };

    self.connectFakeNodes = function() {
        var linkLength = 18;

        // We want to be able to connect all of the fake nodes
        // and create a structure consisting of just them
        var filterOutNonFakeNodes = function(d) {
            return d.nodeType == 'middle';
        }

        var nucsToNodes = {};
        var fakeNodes = self.nodes.filter(filterOutNonFakeNodes);
        var linked = new Set();

        // initialize the nucleotides to nodes
        for (var i = 1; i <= self.nodes.length; i++) 
            nucsToNodes[i] = [];

        for (i = 0; i < fakeNodes.length; i++) {
            var thisNode = fakeNodes[i];

            // each fake node represents a certain set of nucleotdies (thisNode.nucs)
            for (var j = 0; j < thisNode.nucs.length; j++) {
                var thisNuc = thisNode.nucs[j];

                // check to see if this nucleotide has been seen in another fake node
                // if it has, then we add a link between the two nodes
                for (var k = 0; k < nucsToNodes[thisNuc].length; k++) {
                    if (linked.has(JSON.stringify([nucsToNodes[thisNuc][k].uid, thisNode.uid].sort())))
                        continue; //already linked

                    var distance = nucsToNodes[thisNuc][k].radius + thisNode.radius;

                    self.links.push({"source": nucsToNodes[thisNuc][k],
                                      "target": thisNode,
                                      "value": distance / linkLength,
                                      "linkType": "fake_fake"});

                    // note that we've already seen this link
                    linked.add(JSON.stringify([nucsToNodes[thisNuc][k].uid, thisNode.uid].sort()));
                }

                nucsToNodes[thisNuc].push(thisNode);
            }
        }

        return self;

    };

    self.elementsToJson = function() {
        /* Convert a set of secondary structure elements to a json
         * representation of the graph that can be used with d3's
         * force-directed layout to generate a visualization of 
         * the structure.
         */
        pt = self.pairtable;
        elements = self.elements;

        self.nodes = [];
        self.links = [];

        //create a reverse lookup so we can find out the type
        //of element that a node is part of
        elemTypes = {};

        //sort so that we count stems last
        self.elements.sort();

        for (var i = 0; i < self.elements.length; i++) {
            nucs = self.elements[i][2];
            for (j = 0; j < nucs.length; j++) {
                elemTypes[nucs[j]] = self.elements[i][0];
            }
        }

        for (i = 1; i <= pt[0]; i++) {
            //create a node for each nucleotide
            self.nodes.push({'name': self.seq[i-1],
                             'num': i,
                             'radius': 6,
                             'rna': self,
                             'nodeType': 'nucleotide',
                             'structName': self.structName,
                             'elemType': elemTypes[i],
                             'uid': generateUUID() });
        }


        for (i = 1; i <= pt[0]; i++) {

            if (pt[i] !== 0) {
                // base-pair links
                self.links.push({'source': self.nodes[i-1],
                                 'target': self.nodes[pt[i]-1],
                                 'linkType': 'basepair',
                                 'value': 1,
                                 'uid': generateUUID() });
            }

            if (i > 1) {
                // backbone links
                self.links.push({'source': self.nodes[i-2],
                                 'target': self.nodes[i-1],
                                 'linkType': 'backbone',
                                 'value': 1,
                                 'uid': generateUUID() });
            }
        }

        //add the pseudoknot links
        for (i = 0; i < self.pseudoknotPairs.length; i++) {
            self.links.push({'source': self.nodes[self.pseudoknotPairs[i][0]-1],
                            'target': self.nodes[self.pseudoknotPairs[i][1]-1],
                            'linkType': 'pseudoknot',
                            'value': 1,
                            'uid': generateUUID()});
        }

        if (self.circular) {
            self.links.push({'source': self.nodes[0],
                            'target': self.nodes[self.rnaLength-1],
                            'linkType': 'backbone',
                            'value': 1,
                            'uid': generateUUID() });

        }

        return self;
    };

    self.ptToElements = function(pt, level, i, j) {
        /* Convert a pair table to a list of secondary structure 
         * elements:
         *
         * [['s',1,[2,3]]
         *
         * The 's' indicates that an element can be a stem. It can also be
         * an interior loop ('i'), a hairpin loop ('h') or a multiloop ('m')
         *
         * The second number (1 in this case) indicates the depth or
         * how many base pairs have to be broken to get to this element.
         *
         * Finally, there is the list of nucleotides which are part of
         * of this element.
         */
        var elements = [];
        var u5 = [i-1];
        var u3 = [j+1];

        if (i > j)
            return [];
            
            //iterate over the unpaired regions on either side
            //this is either 5' and 3' unpaired if level == 0
            //or an interior loop or a multiloop
            for (; pt[i] === 0; i++) { u5.push(i); }
            for (; pt[j] === 0; j--) { u3.push(j); }

            if (i > j) {
                //hairpin loop or one large unpaired molecule
                u5.push(i);
                if (level === 0)
                    return [['e',level, u5.sort(numberSort)]];
                else
                    return [['h',level, u5.sort(numberSort)]];
            }

            if (pt[i] != j) {
                //multiloop
                var m = u5;
                var k = i;

                // the nucleotide before and the starting nucleotide
                m.push(k);
                while (k <= j) {
                    // recurse into a stem
                    elements = elements.concat(self.ptToElements(pt, level, k, pt[k]));

                    // add the nucleotides between stems
                    m.push(pt[k]);
                    k = pt[k] + 1;
                    for (; pt[k] === 0 && k <= j; k++) { m.push(k);}
                    m.push(k);
                }
                m.pop();
                m = m.concat(u3);
                
                if (m.length > 0) {
                    if (level === 0)
                        elements.push(['e', level, m.sort(numberSort)]);
                    else
                        elements.push(['m', level, m.sort(numberSort)]);
                }
                
                return elements;
            }

            if (pt[i] === j) {
                //interior loop
                u5.push(i);
                u3.push(j);

                combined = u5.concat(u3);
                if (combined.length > 4) {
                    if (level === 0)
                        elements.push(['e',level, u5.concat(u3).sort(numberSort)]);
                    else
                        elements.push(['i',level, u5.concat(u3).sort(numberSort)]);
                }
            } 

            var s = [];
            //go through the stem
            while (pt[i] === j && i < j) {
                //one stem
                s.push(i);
                s.push(j);

                i += 1;
                j -= 1;

                level += 1;
            }

            u5 = [i-1];
            u3 = [j+1];
            elements.push(['s', level, s.sort(numberSort)]);

        return elements.concat(self.ptToElements(pt, level, i, j));
    };

    self.addLabels = function(startNumber, labelInterval) {
        if (arguments.length  === 0) {
            startNumber = 1;
            labelInterval = 10;
        }

        if (arguments.length === 1) 
            labelInterval = 10;

        if (labelInterval === 0)
            return self;

        if (labelInterval <= 0) 
            console.log('The label interval entered in invalid:', labelInterval);

        for (i = 1; i <= pt[0]; i++) {
            // add labels
            if (i % labelInterval === 0) {
                //create a node for each label
                var newX, newY;

                thisNode = self.nodes[i-1]

                if (self.rnaLength == 1) {
                    nextVec = [thisNode.x - 15, thisNode.y]
                    prevVec = [thisNode.x - 15, thisNode.y]
                } else {
                    // if we're labelling the first node, then label it in relation to the last
                    if (i == 1)
                        prevNode = self.nodes[self.rnaLength - 1];
                    else
                        prevNode = self.nodes[i - 2];

                    // if we're labelling the last node, then label it in relation to the first
                    if (i == self.rnaLength)
                        nextNode = self.nodes[0];
                    else
                        nextNode = self.nodes[i];

                    // this nucleotide and its neighbors are paired
                    if (self.pairtable[nextNode.num] !== 0 &&
                        self.pairtable[prevNode.num] !== 0 &&
                        self.pairtable[thisNode.num] !== 0) {
                        prevNode = nextNode = self.nodes[self.pairtable[thisNode.num]-1]
                    }

                    // this node is paired but at least one of its neighbors is unpaired
                    // place the label in the direction of the two neighbors
                    if (self.pairtable[thisNode.num] !== 0 && (
                        self.pairtable[nextNode.num] === 0 ||
                        self.pairtable[prevNode.num] === 0)) {
                        nextVec = [thisNode.x - nextNode.x, thisNode.y - nextNode.y];
                        prevVec = [thisNode.x - prevNode.x, thisNode.y - prevNode.y];

                    } else {
                        nextVec = [nextNode.x - thisNode.x, nextNode.y - thisNode.y];
                        prevVec = [prevNode.x - thisNode.x, prevNode.y - thisNode.y];
                    }
                }

                combinedVec = [nextVec[0] + prevVec[0], nextVec[1] + prevVec[1]];
                vecLength = Math.sqrt(combinedVec[0] * combinedVec[0] + combinedVec[1] * combinedVec[1]);
                normedVec = [combinedVec[0] / vecLength, combinedVec[1] / vecLength];
                offsetVec = [-15 * normedVec[0], -15 * normedVec[1]];

                newX = self.nodes[i-1].x + offsetVec[0];
                newY = self.nodes[i-1].y + offsetVec[1];

                newNode = {'name': i + startNumber - 1,
                                 'num': -1,
                                 'radius': 6,
                                 'rna': self,
                                 'nodeType': 'label',
                                 'structName': self.structName,
                                 'elemType': 'l',
                                 'x': newX,
                                 'y': newY,
                                 'px': newX,
                                 'py': newY,
                                 'uid': generateUUID() };
                newLink = {'source': self.nodes[i-1],
                            'target': newNode,
                            'value': 1,
                            'linkType': 'label_link',
                            'uid': generateUUID() };

                self.nodes.push(newNode);
                self.links.push(newLink);
            }
        }

        return self;
    };

    self.recalculateElements = function() {
        self.removePseudoknots();
        self.elements = self.ptToElements(self.pairtable, 0, 1, self.dotbracket.length);

        if (self.circular) {
            //check to see if the external loop is a hairpin or a multiloop
            externalLoop = self.elements.filter(function(d) { if (d[0] == 'e') return true; });

            if (externalLoop.length > 0) {
                eloop = externalLoop[0];
                nucs = eloop[2].sort(numberSort);

                prev = nucs[0];
                hloop = true;
                numGreater = 0;
                for (var i = 1; i < nucs.length; i++) {
                    if (nucs[i] - prev > 1) {
                        numGreater += 1;
                    }
                    prev = nucs[i];
                }

                if (numGreater == 1) {
                    eloop[0] = 'h';
                } else if (numGreater == 2) {
                    eloop[0] = 'i';
                } else {
                    eloop[0] = 'm';
                }
            }
        }

        return self;
    };

    self.removePseudoknots = function() {
        if (self.pairtable.length > 1)
            self.pseudoknotPairs = self.pseudoknotPairs.concat(rnaUtilities.removePseudoknotsFromPairtable(self.pairtable));

        return self;
    };

    self.addPseudoknots = function() {
        /* Add all of the pseudoknot pairs which are stored outside
         * of the pairtable back to the pairtable
         */
        var pt = self.pairtable;
        var pseudoknotPairs = self.pseudoknotPairs;

        for (i = 0; i < pseudoknotPairs.length; i++) {
            pt[pseudoknotPairs[i][0]] = pseudoknotPairs[i][1];
            pt[pseudoknotPairs[i][1]] = pseudoknotPairs[i][0];
        }

        self.pseudoknotPairs = [];
        return self;
    };

    if (self.rnaLength > 0)
        self.recalculateElements();
}

moleculesToJson = function(moleculesJson) {
    /* Convert a list of RNA and protein molecules to a list of RNAGraph
     * ProteinGraph and extraLinks structure */

    var nodes = {}; //index the nodes by uid
    var graphs = [];
    var extraLinks = [];


    // Create the graphs for each molecule
    for (var i = 0; i < moleculesJson.molecules.length; i++) {
        var molecule = moleculesJson.molecules[i];

        if (molecule.type == 'rna') {
            rg = new RNAGraph(molecule.seq, molecule.ss, molecule.header);
            rg.circularizeExternal = true;
            rg.elementsToJson()
            .addPositions('nucleotide', molecule.positions)
            .addLabels()
            .reinforceStems()
            .reinforceLoops();

            
        } else if (molecule.type == 'protein') {
            rg = new ProteinGraph(molecule.header, molecule.size);

        }

        rg.addUids(molecule.uids);

        for (var j = 0; j < rg.nodes.length; j++) {
            nodes[rg.nodes[j].uid] = rg.nodes[j];
        }

        graphs.push(rg);
    }

    //Add the extra links
    for (i = 0; i < moleculesJson.extraLinks.length; i++) {
        link = moleculesJson.extraLinks[i];
        
        link.source = nodes[link.source];
        link.target = nodes[link.target];
        link.uid = generateUUID();

        extraLinks.push(link);
    }

    return {"graphs": graphs, "extraLinks": extraLinks};
}
var number_sort = function(a,b) { return a - b; };

function RNAUtilities() {
    var self = this;

    // the brackets to use when constructing dotbracket strings
    // with pseudoknots
    self.bracket_left =  "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("");
    self.bracket_right = ")]}>abcdefghijklmnopqrstuvwxyz".split("");

    self.inverse_brackets = function(bracket) {
        res = {};
        for (i = 0; i < bracket.length; i++) {
            res[bracket[i]] = i;
        }
        return res;
    };

    self.maximumMatching = function maximumMatching(pt){
        // Courtesy of the great Ronny Lorenz

        var n = pt[0];
        var TURN = 0;    //minimal number of nucleotides in the hairpin

        /* array init */
        mm = new Array(n + 1);
        for(var i = 0; i <= n; i++){
            mm[i] = new Array(n + 1);
            for(var j = i; j <= n; j++)
            mm[i][j] = 0;
        }
        var maximum = 0;

        /* actual computation */
        for(var i = n - TURN - 1; i > 0; i--)

        for(var j = i + TURN + 1; j <= n; j++){
            maximum = mm[i][j-1];

            for(var l = j - TURN - 1; l >= i; l--) {
                if(pt[l] === j) {

                    // we have a base pair here
                    maximum = Math.max(maximum, ((l > i) ? mm[i][l-1] : 0) + 1 + ((j - l - 1 > 0) ? mm[l+1][j-1] : 0));
                }
            }

            mm[i][j] = maximum;
        }

        maximum = mm[1][n];

        return mm;
    };

    self.backtrackMaximumMatching = function(mm, old_pt) {
      var pt = Array.apply(null, 
                           Array(mm.length)).map(function() { return 0 }); 
                           //create an array containing zeros

      self.mm_bt(mm, pt, old_pt, 1, mm.length-1);
      return pt;
    }

    self.mm_bt = function(mm, pt, old_pt, i, j){
        // Create a pairtable from the backtracking
      var maximum = mm[i][j];
      var TURN = 0;

      if(j - i - 1 < TURN) return;    /* no more pairs */

      if(mm[i][j-1] == maximum){      /* j is unpaired */
        self.mm_bt(mm, pt, old_pt, i, j-1);
        return;
      }

      for(var q = j - TURN - 1; q >= i; q--){  /* j is paired with some q */
        if (old_pt[j] !== q)
            continue;

        var left_part     = (q > i) ? mm[i][q-1] : 0;
        var enclosed_part = (j - q - 1 > 0) ? mm[q+1][j-1] : 0;

        if(left_part + enclosed_part + 1 == maximum) {
            // there's a base pair between j and q
            pt[q] = j;
            pt[j] = q;

            if(i < q) 
                self.mm_bt(mm, pt, old_pt, i, q - 1);

            self.mm_bt(mm, pt, old_pt, q + 1, j - 1);
            return;
        }
      }

      //alert(i + "," + j + ": backtracking failed!");
      console.log("FAILED!!!" + i + "," + j + ": backtracking failed!");

    };

    self.dotbracketToPairtable = function(dotbracket) {
        // create an array and initialize it to 0
        pt = Array.apply(null, new Array(dotbracket.length + 1)).map(Number.prototype.valueOf,0);
        
        //  the first element is always the length of the RNA molecule
        pt[0] = dotbracket.length;

        // store the pairing partners for each symbol
        stack = {};
        for (i = 0; i < self.bracket_left.length; i++) {
            stack[i] = [];
        }

        // lookup the index of each symbol in the bracket array
        inverse_bracket_left = self.inverse_brackets(self.bracket_left);
        inverse_bracket_right = self.inverse_brackets(self.bracket_right);

        for (i = 0; i < dotbracket.length; i++) {
            a = dotbracket[i];
            ni = i + 1;

            if (a == '.') {
                // unpaired
                pt[ni] = 0;
            } else {
                if (a in inverse_bracket_left) {
                    // open pair?
                    stack[inverse_bracket_left[a]].push(ni);
                } else if (a in inverse_bracket_right){
                    // close pair?
                    j = stack[inverse_bracket_right[a]].pop();

                    pt[ni] = j;
                    pt[j] = ni;
                } else {
                    throw "Unknown symbol in dotbracket string";
                }
            }
        }

        for (key in stack) {
            if (stack[key].length > 0) {
                throw "Unmatched base at position " + stack[key][0];
            }
        }

        return pt;
    };

    self.insert_into_stack = function(stack, i, j) {
        var k = 0;
        while (stack[k].length > 0 && stack[k][stack[k].length - 1] < j) {
            k += 1;
        }

        stack[k].push(j);
        return k;
    };

    self.delete_from_stack = function(stack, j) {
        var k = 0;
        while (stack[k].length === 0 || stack[k][stack[k].length-1] != j) {
            k += 1;
        }
        stack[k].pop();
        return k;
    };

    self.pairtableToDotbracket = function(pt) {
        // store the pairing partners for each symbol
        stack = {};
        for (i = 0; i < pt[0]; i++) {
            stack[i] = [];
        }

        seen = {};
        res = "";
        for (i = 1; i < pt[0] + 1; i++) {
            if (pt[i] !== 0 && pt[i] in seen) {
                throw "Invalid pairtable contains duplicate entries";
            }
            seen[pt[i]] = true;

            if (pt[i] === 0) {
                res += '.';
            } else {
                if (pt[i] > i) {
                    res += self.bracket_left[self.insert_into_stack(stack, i, pt[i])];
                } else {
                    res += self.bracket_right[self.delete_from_stack(stack, i)];
                }
            }
        }

        return res;
    };

    self.find_unmatched = function(pt, from, to) {
        /*
         * Find unmatched nucleotides in this molecule.
         */
        var to_remove = [];
        var unmatched = [];

        var orig_from = from;
        var orig_to = to;

        for (var i = from; i <= to; i++)
            if (pt[i] !== 0 && (pt[i] < from || pt[i] > to))
                unmatched.push([i,pt[i]]);

        for (i = orig_from; i <= orig_to; i++) {
            while (pt[i] === 0 && i <= orig_to) i++;

            to = pt[i];

            while (pt[i] === to) {
                i++;
                to--;
            }
            
            to_remove = to_remove.concat(self.find_unmatched(pt, i, to));
        }

        if (unmatched.length > 0)
            to_remove.push(unmatched);

        return to_remove;
    };

    self.removePseudoknotsFromPairtable = function(pt) {
        /* Remove the pseudoknots from this structure in such a fashion
         * that the least amount of base-pairs need to be broken
         *
         * The pairtable is manipulated in place and a list of tuples
         * indicating the broken base pairs is returned.
         */

        var mm = self.maximumMatching(pt);
        var new_pt = self.backtrackMaximumMatching(mm, pt);
        var removed = [];

        for (var i = 1; i < pt.length; i++) {
            if (pt[i] < i)
                continue;

            if (new_pt[i] != pt[i])  {
                removed.push([i, pt[i]]);
                pt[pt[i]] = 0;
                pt[i] = 0;
            }
        }

        return removed;
    };

}

rnaUtilities = new RNAUtilities();
simpleXyCoordinates = function(pair_table)
{
  var INIT_ANGLE=0.;     /* initial bending angle */
  var INIT_X = 100.;     /* coordinate of first digit */
  var INIT_Y = 100.;     /* see above */
  var RADIUS =  15.;

  var x = [], y = [];

  var i, len;
  var  alpha;

  len = pair_table[0];
  var angle = Array.apply(null, new Array(len+5)).map(Number.prototype.valueOf,0); 
  var loop_size = Array.apply(null, new Array(16+Math.floor(len/5)))
                    .map(Number.prototype.valueOf, 0); 
  var stack_size = Array.apply(null, new Array(16+Math.floor(len/5)))
                    .map(Number.prototype.valueOf, 0); 

  lp = stk = 0;
  var PIHALF = Math.PI / 2;


  loop = function(i, j, pair_table)
  /* i, j are the positions AFTER the last pair of a stack; i.e
     i-1 and j+1 are paired. */
  {
      var count = 2;   /* counts the VERTICES of a loop polygon; that's
                          NOT necessarily the number of unpaired bases!
                          Upon entry the loop has already 2 vertices, namely
                          the pair i-1/j+1.  */

  var    r = 0, bubble = 0; /* bubble counts the unpaired digits in loops */

  var    i_old, partner, k, l, start_k, start_l, fill, ladder;
  var    begin, v, diff;
  var  polygon;

  var remember = Array.apply(null, new Array((1+Math.floor((j-i)/5)*2))).map(Number.prototype.valueOf, 0);

  i_old = i-1, j++;         /* j has now been set to the partner of the
                               previous pair for correct while-loop
                               termination.  */
  while (i != j) {
      partner = pair_table[i];
      if ((!partner) || (i==0))
          i++, count++, bubble++;
      else {
          count += 2;
          k = i, l = partner;    /* beginning of stack */
          remember[++r] = k;
          remember[++r] = l;
          i = partner+1;         /* next i for the current loop */

          start_k = k, start_l = l;
          ladder = 0;
          do {
              k++, l--, ladder++;        /* go along the stack region */
          }
          while (pair_table[k] == l);

          fill = ladder-2;
          if (ladder >= 2) {
              angle[start_k+1+fill] += PIHALF;   /*  Loop entries and    */
              angle[start_l-1-fill] += PIHALF;   /*  exits get an        */
              angle[start_k]        += PIHALF;   /*  additional PI/2.    */
              angle[start_l]        += PIHALF;   /*  Why ? (exercise)    */
              if (ladder > 2) {
                  for (; fill >= 1; fill--) {
                      angle[start_k+fill] = Math.PI;    /*  fill in the angles  */
                      angle[start_l-fill] = Math.PI;    /*  for the backbone    */
                  }
              }
          }
          stack_size[++stk] = ladder;
          loop(k, l, pair_table);
      }
  }

  polygon = Math.PI*(count-2)/count; /* bending angle in loop polygon */
  remember[++r] = j;
  begin = i_old < 0 ? 0 : i_old;
  for (v = 1; v <= r; v++) {
      diff  = remember[v]-begin;
      for (fill = 0; fill <= diff; fill++)
      angle[begin+fill] += polygon;
      if (v > r)
          break;
      begin = remember[++v];
  }
  loop_size[++lp] = bubble;
  }

  loop(0, len+1, pair_table);
  loop_size[lp] -= 2;     /* correct for cheating with function loop */

  alpha = INIT_ANGLE;
  x[0]  = INIT_X;
  y[0]  = INIT_Y;

  poss = [];

  poss.push([x[0], y[0]]);
  for (i = 1; i < len; i++) {
      x[i] = x[i-1]+RADIUS*Math.cos(alpha);
      y[i] = y[i-1]+RADIUS*Math.sin(alpha);

      poss.push([x[i], y[i]]);
      alpha += Math.PI-angle[i+1];
  }

  return poss;
}
