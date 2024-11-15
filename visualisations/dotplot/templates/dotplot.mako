<!DOCTYPE HTML>
<%
    hdadict = trans.security.encode_dict_ids( hda.to_dict() )
%>
<html>
<head>
<title>${hda.name} | ${visualization_name}</title>

<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<script type="text/javascript" src="../../../../../static/scripts/libs/d3.js"></script>
<script type="text/javascript" src="../../../../../static/scripts/libs/jquery/jquery.js"></script>
<script type="text/javascript" src="../../../../../static/scripts/libs/jquery/jquery.migrate.js"></script>


</head>


<style>

.chart div {
    font: 10px sans-serif;
    background-color: steelblue;
    text-align: right;
    padding: 3px;
    margin: 1px;
    color: white;
}
.chart2 div {
    font: 10px sans-serif;
    background-color: red;
    text-align: right;
    padding: 3px;
    margin: 1px;
    color: white;
}

text.active {
    fill: red;
    font-weight: bolder;
    font-size: 18px;
}

rect.active{
    fill: red;
}

.background {
      fill: #fff;
}

line {
      stroke: #eee;
}

svg {
      font: 10px sans-serif;
}

</style>

<body>
<h1>DotPlot Matrix<h1>

<script type="text/javascript">

function parse_ps(img) {
    var maindic={},
        seqfound=false,
        sequence="",
        bpp=[],
        mlp=[];
    var keys=["source","target","value"];

    var dic= img.split("\n");
    for (var i=0;i<dic.length ; i++){
        f=dic[i];
        if (seqfound){
            sequence=f.substring(0,f.length-1);
            seqfound=false;
        }
        if (f.search("/sequence") != -1){
            seqfound=true;
            //console.log(f);
        }
        if (f.search(" ubox") != -1 && i>15){
            var a=f.split(" ")[0];
            b=f.split(" ")[1];
            c=f.split(" ")[2];
            var row = {};
            row[keys[0]]=a;
            row[keys[1]]=b;
            row[keys[2]]=c;
            bpp.push(row);
        }
        if (f.search(" lbox") != -1 && i>15){
            a=f.split(" ")[0];
            b=f.split(" ")[1];
            var row = {};
            row[keys[0]]=a;
            row[keys[1]]=b;
            mlp.push(row);
        }
    }
    maindic["base-pairing-probabilities"]=bpp;
    maindic["optimal-structure"]=mlp;
    maindic["sequence"]=sequence;
    
    return JSON.stringify(maindic);
}

function mouseover(p) {
    d3.selectAll(".rowtext").classed("active", function(d, i) { return i == p.y; });
    d3.selectAll(".coltext").classed("active", function(d, i) { return i == p.x; });
    d3.selectAll(".row .cell rect").classed("active", function(d, i) { return d.y == p.y && d.x == p.x; });
    d3.selectAll(".row .cell text").attr("display",function(d, i) { if (d.y == p.y && d.x == p.x) return "true"; else return "none" });
}

function mouseout() {
    d3.selectAll("text").classed("active", false);
    d3.selectAll(".cell").classed("active", false);
    d3.selectAll(".row .cell text").attr("display","none");
}

function showtext() {
    d3.selectAll(".row .cell text")
    .attr("display","true")
}

function hidetext() {
    d3.selectAll(".row .cell text")
    .attr("display","none")
}


$(document).ready(function() {
    var hdaId   = '${trans.security.encode_id( hda.id )}',
        hdaExt  = '${hda.ext}',
        ajaxUrl = "${h.url_for( controller='/datasets', action='index')}/" + hdaId + "/display?to_ext=" + hdaExt;


    %if hda.ext == 'rna_eps':
    $.ajax( ajaxUrl, {
        success: function(psImg) {
            var psjson=parse_ps(psImg);

        var nodes = [{x: 30, y: 50}];
        var bpm=JSON.parse(psjson);
    
        var margin = {top: 80, right: 80, bottom: 10, left: 80},
            width = 900,
            height = 900;
        var x = d3.scale.ordinal().rangeBands([0, width]),
            z = d3.scale.linear().domain([0, 1]).clamp(true),
            c = d3.scale.category10().domain(d3.range(10));

        var force = d3.layout.force()
            .charge(-120)
            .linkDistance(30)
            .size([width, height]);
    
        var svg = d3.select("body").append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom)
            .style("margin-left", -margin.left/2 + "px")
            .append("g")
            .attr("fill", "black")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
    
        var matrix = [];   
    
        var seq = [bpm.sequence],  
            isequence = bpm.sequence.split(""),
            n = isequence.length;
    
        var seq = bpm["sequence"];
        var isequence = bpm["sequence"];
        var n = seq.length;
        isequence=seq.split("");
    
        x.domain(d3.range(n)); 
        isequence.forEach(function(base, i) {
            base.index = i;
            matrix[i] = d3.range(n).map(function(j) { return {x: j, y: i, z: 0}; });
        });
        bpm["base-pairing-probabilities"].forEach(function(link) {
            matrix[link.source-1][link.target-1].z += link.value;
            //matrix[link.target-1][link.source-1].z += link.value;
        });
    
        bpm["optimal-structure"].forEach(function(link) {
            matrix[link.target-1][link.source-1].z = 1;
        });
    
    
        var column = svg.selectAll(".column")
            .data(matrix)
            .enter().append("g")
            .attr("class", "column")
            .attr("transform", function(d, i) { return "translate(" + x(i) + ")rotate(-90)"; });
    
        column.append("line")
            .attr("x1", -width);
    
        column.append("text")
            .attr("class","coltext")
            .attr("x", 6)
            .attr("y", x.rangeBand() / 2)
            .attr("dy", ".32em")
            .attr("text-anchor", "start")
            .text(function(d, i) { return isequence[i]; });

        var row = svg.selectAll(".row")
            .data(matrix)
            .enter().append("g")
            .attr("class", "row")
            .attr("transform", function(d, i) { return "translate(0," + x(i) + ")"; })
            .each(row);
    
        row.append("line")
            .attr("x2", width);
    
        row.append("text")
            .attr("class","rowtext")
            .attr("x", -6)
            .attr("y", x.rangeBand() / 2)
            .attr("dy", ".32em")
            .attr("text-anchor", "end")
            .text(function(d, i) { return isequence[i]; });
    
        function row(row) {
            var cell = d3.select(this).selectAll(".cell")
                .data(row.filter(function(d) { return d.z; }))
                .enter().append("g")
                .attr("class", "cell");
            cell.append("rect")
                .attr("x", function(d) { return x(d.x); })
                .attr("width", x.rangeBand())
                .attr("height", x.rangeBand())
                .style("fill-opacity", function(d) { return z(d.z); })
                .on("mouseover", mouseover)
                .on("mouseout", mouseout);
            cell.append("text")
                .attr("x", function(d) { return x(d.x)+15; })
                .attr("y", x.rangeBand()/2)
                .attr("dy", ".32em")
                .attr("text-anchor", "start")
                .attr("display", "none")
                //.text("hello");
                .text(function(d) { return d.z; })
        }
        }
	});
    %endif
});

</script>


  
</body>
</html>

