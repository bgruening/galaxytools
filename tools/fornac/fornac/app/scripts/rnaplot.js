import {simpleXyCoordinates} from './simplernaplot.js';
import {ProteinGraph, RNAGraph,moleculesToJson} from './rnagraph.js';
import {rnaUtilities,ColorScheme} from 'rnautils';

import '../styles/rnaplot.css';

function isNormalInteger(str) {
    //http://stackoverflow.com/a/10834843/899470
    return /^\+?(0|[1-9]\d*)$/.test(str);
}

if(typeof(String.prototype.trim) === 'undefined')
    {
        String.prototype.trim = function() 
        {
            return String(this).replace(/^\s+|\s+$/g, '');
        };
    }


export function rnaPlot() {
    var options = {
        'width': 400,
        'height': 400,
        'nucleotideRadius': 5,
        'rnaEdgePadding': 0,     // how far the leftmost, rightmost, topmost and bottomost
                                // nucleotides are from the edge of the plot
        'labelInterval': 0,
        'showNucleotideLabels': true,
        'startNucleotideNumber': 1,
        'bundleExternalLinks': false
    };

    var xScale, yScale;

    function createTransformToFillViewport(xValues, yValues, molName='') {
        // create transform that will scale the x and y values so that
        // they fill the available viewport
    
        // find out leftmost, rightmost, topmost, bottommost positions of each
        // nucleotide so that we can create a scale
        var xExtent = d3.extent(xValues);
        var yExtent = d3.extent(yValues); 

        var NAME_OFFSET = 30;
        if (molName != '')
            yExtent[1] += NAME_OFFSET;

        // add the radius of the nucleotides
        xExtent[0] -= options.nucleotideRadius + options.rnaEdgePadding;
        yExtent[0] -= options.nucleotideRadius + options.rnaEdgePadding;

        xExtent[1] += options.nucleotideRadius + options.rnaEdgePadding;
        yExtent[1] += options.nucleotideRadius + options.rnaEdgePadding;

        // find out how wide and height the molecule
        var xRange = xExtent[1] - xExtent[0];
        var yRange = yExtent[1] - yExtent[0];

        // how much wider / taller is it than the available viewport
        var xExtra = xRange - options.width;
        var yExtra = yRange - options.height;

        // once we have a scale for one dimension, we can create the scale for the other
        // keeping the same expansion / shrinking ratio
        function createOtherScale(firstScale, newDomain, newRange) {
            var scaleFactor = (firstScale.range()[1] - firstScale.range()[0]) / 
                              (firstScale.domain()[1] - firstScale.domain()[0]);
            var newWidth = (newDomain[1] - newDomain[0]) * scaleFactor
            var newMargin = ((newRange[1] - newRange[0]) - newWidth) / 2;

            return {'scaleFactor': scaleFactor, 
                    'scale': d3.scale.linear()
                                     .domain(newDomain)
                                     .range([newRange[0] + newMargin, newRange[1] - newMargin])};
        }

        var ret;

        if (xExtra > yExtra) {
            // we have to shrink more in the x-dimension than the y
            xScale = d3.scale.linear()
            .domain(xExtent)
            .range([0, options.width])

            ret = createOtherScale(xScale, yExtent, [0, options.height]);
            yScale = ret.scale;
        } else {
            // we have to shrink more in the x-dimension than the y
            yScale = d3.scale.linear()
            .domain(yExtent)
            .range([0, options.height])

            ret = createOtherScale(yScale, xExtent, [0, options.width]);
            xScale = ret.scale;
        }

        var xOffset = xScale.range()[0] - xScale.domain()[0];
        var yOffset = yScale.range()[0] - yScale.domain()[0];

        return 'translate(' + -(xScale.domain()[0] * ret.scaleFactor - xScale.range()[0]) + 
                  ',' + -(yScale.domain()[0] * ret.scaleFactor - yScale.range()[0]) + ')' + 
            'scale(' + ret.scaleFactor + ')';
    }

    function createNucleotides(selection, nucleotideNodes) {
        // create groupings for each nucleotide and label
        var gs = selection
        .selectAll('.rna-base')
        .data(nucleotideNodes)
        .enter()
        .append('svg:g')
        .attr('transform', function(d) { 
            return 'translate(' + d.x + ',' + d.y + ')'; 
        });

        var circles = gs.append('svg:circle')
        .attr('r', options.nucleotideRadius)
        .classed('rna-base', true)

        if (options.showNucleotideLabels) {
            var nucleotideLabels = gs.append('svg:text')
            .text(function(d) { return d.name; })
            .attr('text-anchor', 'middle')
            .attr('dominant-baseline', 'central')
            .classed('nucleotide-label', true)
            .append('svg:title')
            .text(function(d) { return d.struct_name + ':' + d.num; });
        }
    }

    function createLabels(selection, labelNodes) {
        // create groupings for each nucleotide and label

        var gs = selection 
        .selectAll('.rnaLabel')
        .data(labelNodes)
        .enter()
        .append('svg:g')
        .attr('transform', function(d) { 
            return 'translate(' + d.x + ',' + d.y + ')'; 
        });

        var numberLabels = gs.append('svg:text')
        .text(function(d) { return d.name; })
        .attr('text-anchor', 'middle')
        .attr('font-weight', 'bold')
        .attr('dominant-baseline', 'central')
        .classed('number-label', true);
    }

    function createName(selection, name) {
        selection.append('svg:text')
        .attr('transform', 'translate(' + xScale.invert(options.width / 2) + ',' + yScale.invert(options.height) + ')')
        .attr('dy', -10)
        .classed('rna-name', true)
        .text(name);
    }

    function makeExternalLinksBundle(selection, links) {
        var nodesDict = {};
        var linksList = [];
        links = links.filter(function(d) { return d.linkType == 'correct' || d.linkType == 'incorrect' || d.linkType == 'extra'; });

        selection.selectAll('[link-type=extra]')
        .remove();


        for (var i = 0; i < links.length; i++) {
            if (links[i].source === null || links[i].target === null)
                continue;

            nodesDict[links[i].source.uid] = links[i].source;
            nodesDict[links[i].target.uid] = links[i].target;

            linksList.push({'source': links[i].source.uid, 'target': links[i].target.uid, 'linkType': links[i].linkType, 'extraLinkType': links[i].extraLinkType}) ;
        }

        var fbundling = d3.ForceEdgeBundling().nodes(nodesDict).edges(linksList)
        .compatibility_threshold(0.8).step_size(0.2);
        var results   = fbundling();

        var d3line = d3.svg.line()
            .x(function(d){return d.x;})
            .y(function(d){return d.y;})
            .interpolate('linear');

        for (var i = 0; i < results.length; i++) {
            var edge_subpoint_data = results[i];
            // for each of the arrays in the results 
            // draw a line between the subdivions points for that edge

            selection.append('path').attr('d', d3line(edge_subpoint_data))
            .style('fill', 'none')
            .attr('link-type', function(d) { return linksList[i].linkType; })
            .attr('extra-link-type', function(d) { return linksList[i].extraLinkType; })
            .style('stroke-opacity',0.4); //use opacity as blending
        }
        
    }

    function createLinks(selection, links) {
        links = links.filter(function(d) { return d.source !== null && d.target !== null; });
        var gs = selection.selectAll('.rna-link')
        .data(links)
        .enter()
        .append('svg:line')
        .attr('x1', function(d) { return d.source.x; })
        .attr('x2', function(d) { return d.target.x; })
        .attr('y1', function(d) { return d.source.y; })
        .attr('y2', function(d) { return d.target.y; })
        .attr('link-type', function(d) { return d.linkType; })
        .attr('extra-link-type', function(d) { return d.extraLinkType; })
        .classed('rna-link', true);
    }

    function chart(selection) {
        selection.each(function(data) {
            // data should be a dictionary containing at least a structure
            // and possibly a sequence
            let rg = new RNAGraph(data.sequence, data.structure, data.name)
                    .recalculateElements()
                    .elementsToJson()
                    .addName(data.name);

            data.rnaGraph = rg;
            // calculate the position of each nucleotide
            // the positions of the labels will be calculated in
            // the addLabels function
            var positions = simpleXyCoordinates(rg.pairtable);
            rg.addPositions('nucleotide', positions)
            .reinforceStems()
            .reinforceLoops()
            .addExtraLinks(data.extraLinks)
            .addLabels(options.startNucleotideNumber, options.labelInterval);

            // create a transform that will fit the molecule to the
            // size of the viewport (canvas, svg, whatever)
            var fillViewportTransform = createTransformToFillViewport(
                rg.nodes.map(function(d) { return d.x; }),
                rg.nodes.map(function(d) { return d.y; }));

            var gTransform = d3.select(this)
            .append('g')
            .attr('transform', fillViewportTransform);

            var nucleotideNodes = rg.nodes.filter(function(d) { 
                return d.nodeType == 'nucleotide'; 
            });

            var labelNodes = rg.nodes.filter(function(d) {
                return d.nodeType == 'label';
            });

            var links = rg.links;

            createLinks(gTransform, links);
            createNucleotides(gTransform, nucleotideNodes);            
            createLabels(gTransform, labelNodes);
            createName(gTransform, data.name);

            if (options.bundleExternalLinks) {
                makeExternalLinksBundle(gTransform, links); 
            }

        });
    }

    chart.width = function(_) {
        if (!arguments.length) return options.width;
        options.width = _;
        return chart;
    };

    chart.height = function(_) {
        if (!arguments.length) return options.height;
        options.height = _;
        return chart;
    };

    chart.showNucleotideLabels = function(_) {
        if (!arguments.length) return options.showNucleotideLabels;
        options.showNucleotideLabels = _;
        return chart;
    };

    chart.rnaEdgePadding = function(_) {
        if (!arguments.length) return options.rnaEdgePadding;
        options.rnaEdgePadding = _;
        return chart;
    };

    chart.nucleotideRadius = function(_) {
        if (!arguments.length) return options.nucleotideRadius;
        options.nucleotideRadius = _;
        return chart;
    };

    chart.labelInterval = function(_) {
        if (!arguments.length) return options.labelInterval;
        options.labelInterval = _;
        return chart;
    };

    chart.showNucleotideLabels = function(_) {
        if (!arguments.length) return options.showNucleotideLabels;
        options.showNucleotideLabels = _;
        return chart;
    };

    chart.startNucleotideNumber = function(_) {
        if (!arguments.length) return options.startNucleotideNumber;
        options.startNucleotideNumber = _;
        return chart;
    };

    chart.bundleExternalLinks = function(_) {
        if (!arguments.length) return options.bundleExternalLinks;
        options.bundleExternalLinks = _;
        return chart;
    };

    return chart;
}
var number_sort = function(a,b) { return a - b; };

function RNAUtilities() {
    var self = this;

    // the brackets to use when constructing dotbracket strings
    // with pseudoknots
    self.bracket_left =  '([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ'.split('');
    self.bracket_right = ')]}>abcdefghijklmnopqrstuvwxyz'.split('');

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

      //alert(i + ',' + j + ': backtracking failed!');
      console.log('FAILED!!!' + i + ',' + j + ': backtracking failed!');

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
                    throw 'Unknown symbol in dotbracket string';
                }
            }
        }

        for (key in stack) {
            if (stack[key].length > 0) {
                throw 'Unmatched base at position ' + stack[key][0];
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
        res = '';
        for (i = 1; i < pt[0] + 1; i++) {
            if (pt[i] !== 0 && pt[i] in seen) {
                throw 'Invalid pairtable contains duplicate entries';
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
