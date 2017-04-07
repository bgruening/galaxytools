import {arraysEqual,RNAUtilities,rnaUtilities} from 'rnautils';

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

if(typeof(String.prototype.trim) === 'undefined')
    {
        String.prototype.trim = function() 
        {
            return String(this).replace(/^\s+|\s+$/g, '');
        };
    }


export function ProteinGraph(structName, size, uid) {
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

export function RNAGraph(seq, dotbracket, structName, startNumber) {
    var self = this;

    self.type = 'rna';
    self.circularizeExternal = false;

    if (arguments.length === 0) {
        self.seq = '';
        self.dotbracket = '';
        self.structName = '';
    } else {
        self.seq = seq;
        self.dotbracket = dotbracket;  //i.e. ..((..))..
        self.structName = structName;
    }

    if (arguments.length < 4) {
        startNumber = 1;
    }

    self.circular = false;

    if (self.dotbracket.length > 0 && self.dotbracket[self.dotbracket.length-1] == '*') {
        //circular RNA
        self.dotbracket = self.dotbracket.slice(0, self.dotbracket.length-1);
        self.circular = true;
    }

    self.uid = generateUUID();

    self.elements = [];            //store the elements and the 
                                   //nucleotides they contain
    self.pseudoknotPairs = [];
    self.nucsToNodes = {};

    self.addUids = function(uids) {
        var nucleotideNodes = self.nodes.filter(function(d) { return d.nodeType == 'nucleotide'; });

        for (var i = 0; i < uids.length && i < nucleotideNodes.length; i++)
            nucleotideNodes[i].uid = uids[i];

        return self;
    };

    self.computePairtable = function() {
        self.pairtable = rnaUtilities.dotbracketToPairtable(self.dotbracket);
    };

    self.removeBreaks = function(targetString) {
        // Remove all chain breaks (denoted with a '&', which indicate
        // that the input represents more than one strand)
        var breaks = [];
        var breakIndex = -1;

        while ((breakIndex = targetString.indexOf('&')) >= 0) {
            breaks.push(breakIndex);
            targetString = targetString.substring(0, breakIndex) + 'oo' + targetString.substring(breakIndex+1, targetString.length);

        }

        return {targetString: targetString,  breaks: breaks};
    };

    var ret = self.removeBreaks(self.dotbracket);
    self.dotbracket = ret.targetString;
    self.dotBracketBreaks = ret.breaks;

    ret = self.removeBreaks(self.seq);
    self.seq = ret.targetString;
    self.seqBreaks = ret.breaks;

    self.calculateStartNumberArray = function() {
        self.startNumberArray = [];
        var breaks = 0;

        for (var i = 0; i < self.dotbracket.length; i++) {
            self.startNumberArray.push(startNumber); 

            if (self.dotbracket[i] == 'o') {
                startNumber = -i;
            }
        }

    };

    self.calculateStartNumberArray();

    self.rnaLength = self.dotbracket.length;

    if (!arraysEqual(self.dotBracketBreaks, self.seqBreaks)) {
        console.log('WARNING: Sequence and structure breaks not equal');
        console.log('WARNING: Using the breaks in the structure');
    }
    
    self.computePairtable();

    self.addPositions = function(nodeType, positions) {
        var labelNodes = self.nodes.filter(function(d) { return d.nodeType == nodeType; });

        for  (var i = 0; i < labelNodes.length; i++) {
            labelNodes[i].x = positions[i][0];
            labelNodes[i].px = positions[i][0];
            labelNodes[i].y = positions[i][1];
            labelNodes[i].py = positions[i][1];
        }

        return self;
    };

    self.breakNodesToFakeNodes = function() {
        // convert all the nodes following breaks to fake nodes
        var labelNodes = self.nodes.filter(function(d) { return d.nodeType == 'nucleotide'; });

        // if a node was an artifical break node, convert it to a middle
        for (var i = 0; i < labelNodes.length; i++) {
            if (self.dotbracket[i] == 'o')
                labelNodes[i].nodeType = 'middle';
        }

        for (var i = 0; i < self.elements.length; i++) {
            var broken = false;

            // change the elemType of the other nodes in the element containing
            // the break
            for (var j = 0; j < self.elements[i][2].length; j++) {
                if (self.dotBracketBreaks.indexOf(self.elements[i][2][j]) >= 0)
                    broken = true;
            }

            if (broken) {
                self.elements[i][2].map(function(x) {
                    if (x == 0)
                        return;
                    self.nodes[x-1].elemType = 'e';
                });
            } else {
                self.elements[i][2].map(function(x) {
                    if (x == 0)
                        return;
                    self.nodes[x-1].elemType = self.elements[i][0];
                });
            }
        }
        return self;
    };

    self.getPositions = function(nodeType) {
        var positions = [];
        var nucleotideNodes = self.nodes.filter(function(d) { return d.nodeType == nodeType; });

        for (var i = 0; i < nucleotideNodes.length; i++)
            positions.push([nucleotideNodes[i].x, nucleotideNodes[i].y]);

        return positions;
    };

    self.getUids = function() {
        /* Get the positions of each node so that they
         * can be passed to elementsToJson later
         */
        var uids = [];
        for (var i = 0; i < self.dotbracket.length; i++)
            uids.push(self.nodes[i].uid);

        return uids;
    };

    self.reinforceStems = function() {
        var pt = self.pairtable;
        var relevantElements = self.elements.filter( function(d) {
            return d[0] == 's' && d[2].length >= 4;
        });

        for (var i = 0; i < relevantElements.length; i++) {
            var allNucs = relevantElements[i][2];
            var nucs = allNucs.slice(0, allNucs.length / 2);

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

        for (var i = 0; i < self.elements.length; i++) {
            if (self.elements[i][0] == 's' || (!self.circularizeExternal && self.elements[i][0] == 'e'))
                continue;

            var nucs = self.elements[i][2].filter(filterNucs);

            if (self.elements[i][0] == 'e') {
                var newNode1 = {'name': '',
                    'num': -3,
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
                    'num': -2,
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
    };

    self.addFakeNode = function(nucs) {
        var linkLength = 18; //make sure this is consistent with the value in force.js
        var nodeWidth = 6;
        var angle = (3.1415 * 2) / (2 * nucs.length);
        var radius =  linkLength / (2 * Math.tan(angle));

        var fakeNodeUid = '';

        for (var i = 0; i < nucs.length; i++)
            fakeNodeUid += self.nodes[nucs[i]-1].uid;

        var newNode = {'name': '',
                         'num': -1,
                         //'radius': 18 * radius -6,
                         'radius': radius,
                         'rna': self,
                         'nodeType': 'middle',
                         'elemType': 'f',
                         'nucs': nucs,
                         'uid': fakeNodeUid };
        self.nodes.push(newNode);

        var newX = 0;
        var newY = 0;
        var coordsCounted = 0;

        angle = (nucs.length - 2) * 3.14159 / (2 * nucs.length);
        radius = 0.5 / Math.cos(angle);

        for (var j = 0; j < nucs.length; j++) {
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

            var ia = ((nucs.length - 2) * 3.14159) / nucs.length;
            var c = 2 * Math.cos(3.14159 / 2 - ia / 2);
            //link to over-neighbor
            self.links.push({'source': self.nodes[nucs[j] - 1],
                             'target': self.nodes[nucs[(j + 2) % nucs.length] - 1],
                             'linkType': 'fake',
                             'value': c});

            // calculate the mean of the coordinats in this loop
            // and place the fake node there
            var fromNode = self.nodes[nucs[j]-1];
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

        for (var i = 0; i < fakeNodes.length; i++) {
            var thisNode = fakeNodes[i];

            // each fake node represents a certain set of nucleotides (thisNode.nucs)
            for (var j = 0; j < thisNode.nucs.length; j++) {
                var thisNuc = thisNode.nucs[j];

                // check to see if this nucleotide has been seen in another fake node
                // if it has, then we add a link between the two nodes
                for (var k = 0; k < nucsToNodes[thisNuc].length; k++) {
                    if (linked.has(JSON.stringify([nucsToNodes[thisNuc][k].uid, thisNode.uid].sort())))
                        continue; //already linked

                    var distance = nucsToNodes[thisNuc][k].radius + thisNode.radius;

                    self.links.push({'source': nucsToNodes[thisNuc][k],
                                      'target': thisNode,
                                      'value': distance / linkLength,
                                      'linkType': 'fake_fake'});

                    // note that we've already seen this link
                    linked.add(JSON.stringify([nucsToNodes[thisNuc][k].uid, thisNode.uid].sort()));
                }

                nucsToNodes[thisNuc].push(thisNode);
            }
        }

        return self;

    };

    self.addExtraLinks = function(extraLinks) {
        if (typeof extraLinks == 'undefined')
            return self;

        for (var i = 0; i < extraLinks.length; i++) {
            var source = self.getNodeFromNucleotides(extraLinks[i].from);
            var target = self.getNodeFromNucleotides(extraLinks[i].to);

            var newLink = {'source': source, 'target': target, 'linkType': 'extra',
                'extraLinkType': extraLinks[i].linkType, 'uid': generateUUID() };

                self.links.push(newLink);
        }

        return self;
    }


    self.elementsToJson = function() {
        /* Convert a set of secondary structure elements to a json
         * representation of the graph that can be used with d3's
         * force-directed layout to generate a visualization of 
         * the structure.
         */
        var pt = self.pairtable;
        var elements = self.elements;

        self.nodes = [];
        self.links = [];

        //create a reverse lookup so we can find out the type
        //of element that a node is part of
        var elemTypes = {};

        //sort so that we count stems last
        self.elements.sort();

        for (var i = 0; i < self.elements.length; i++) {
            var nucs = self.elements[i][2];
            for (var j = 0; j < nucs.length; j++) {
                elemTypes[nucs[j]] = self.elements[i][0];
            }
        }

        for (var i = 1; i <= pt[0]; i++) {
            var nodeName = self.seq[i-1];

            if (self.dotBracketBreaks.indexOf(i-1) >= 0 ||
                self.dotBracketBreaks.indexOf(i-2) >= 0) {
                nodeName = '';
            }

            //create a node for each nucleotide
            self.nodes.push({'name': nodeName,
                             'num': i + self.startNumberArray[i-1] - 1,
                             'radius': 5,
                             'rna': self,
                             'nodeType': 'nucleotide',
                             'structName': self.structName,
                             'elemType': elemTypes[i],
                             'uid': generateUUID(),
                             'linked': false});
        }

        for (var i = 0; i < self.nodes.length; i++) {
            if (i === 0) 
                self.nodes[i].prevNode = null;
            else {
                self.nodes[i].prevNode = self.nodes[i-1];
            }

            if (i == self.nodes.length-1) 
                self.nodes[i].nextNode = null;
            else {
                self.nodes[i].nextNode = self.nodes[i+1];
            }
        }

        for (var i = 1; i <= pt[0]; i++) {

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
                if (self.dotBracketBreaks.indexOf(i-1) === -1 &&
                    self.dotBracketBreaks.indexOf(i-2) == -1 &&
                    self.dotBracketBreaks.indexOf(i-3) == -1) {
                    // there is no break in the strands here
                    // we can add a backbone link
                    self.links.push({'source': self.nodes[i-2],
                                    'target': self.nodes[i-1],
                                    'linkType': 'backbone',
                                    'value': 1,
                                    'uid': generateUUID() });
                    self.nodes[i-1].linked = true;
                }
            }
        }

        //add the pseudoknot links
        for (var i = 0; i < self.pseudoknotPairs.length; i++) {
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
                else {
                    // check to see if we have chain breaks due
                    // to multiple strands in the input
                    var external = false;
                    var left = [];
                    var right = [];
                    for (var k = 0; k < u5.length; k++) {
                        if (external)
                            right.push(u5[k]);
                        else
                            left.push(u5[k]);

                        if (self.dotBracketBreaks.indexOf(u5[k]) >= 0)
                            external = true;
                    }

                    if (external) {
                        return [['h',level, u5.sort(numberSort)]];
                    }
                    else
                        // if not, this is a simple hairpin loop
                        return [['h',level, u5.sort(numberSort)]];
                }
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

                var combined = u5.concat(u3);
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

        for (var i = 1; i <= self.pairtable[0]; i++) {
            // add labels
            if (i % labelInterval === 0) {
                //create a node for each label
                var newX, newY;

                var thisNode = self.nodes[i-1];
                var prevNode, nextNode;
                var prevVec, nextVec;

                if (self.rnaLength == 1) {
                    nextVec = [thisNode.x - 15, thisNode.y];
                    prevVec = [thisNode.x - 15, thisNode.y];
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
                        prevNode = nextNode = self.nodes[self.pairtable[thisNode.num]-1];
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

                var combinedVec = [nextVec[0] + prevVec[0], nextVec[1] + prevVec[1]];
                var vecLength = Math.sqrt(combinedVec[0] * combinedVec[0] + combinedVec[1] * combinedVec[1]);
                var normedVec = [combinedVec[0] / vecLength, combinedVec[1] / vecLength];
                var offsetVec = [-15 * normedVec[0], -15 * normedVec[1]];

                var newX = self.nodes[i-1].x + offsetVec[0];
                var newY = self.nodes[i-1].y + offsetVec[1];

                var newNode = {'name': i + self.startNumberArray[i-1] - 1,
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
                var newLink = {'source': self.nodes[i-1],
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

    self.reassignLinkUids = function() {
        // reassign uids to the links, corresponding to the uids of the two nodes
        // they connect
        var i;

        for (var i = 0; i < self.links.length; i++) {
            self.links[i].uid = self.links[i].source.uid + self.links[i].target.uid;
        }

        return self;
    }

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

        for (var i = 0; i < pseudoknotPairs.length; i++) {
            pt[pseudoknotPairs[i][0]] = pseudoknotPairs[i][1];
            pt[pseudoknotPairs[i][1]] = pseudoknotPairs[i][0];
        }

        self.pseudoknotPairs = [];
        return self;
    };

    self.addName = function(name) {
        if (typeof name == 'undefined') {
            self.name = '';
            return self;
        } else {
            self.name = name;
            return self;
        }
    };


    if (self.rnaLength > 0)
        self.recalculateElements();
}

export function moleculesToJson(moleculesJson) {
    /* Convert a list of RNA and protein molecules to a list of RNAGraph
     * ProteinGraph and extraLinks structure */

    var nodes = {}; //index the nodes by uid
    var graphs = [];
    var extraLinks = [];


    // Create the graphs for each molecule
    for (var i = 0; i < moleculesJson.molecules.length; i++) {
        var molecule = moleculesJson.molecules[i];
        var rg;

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
    for (var i = 0; i < moleculesJson.extraLinks.length; i++) {
        link = moleculesJson.extraLinks[i];
        
        link.source = nodes[link.source];
        link.target = nodes[link.target];
        link.uid = generateUUID();

        extraLinks.push(link);
    }

    return {'graphs': graphs, 'extraLinks': extraLinks};
};
