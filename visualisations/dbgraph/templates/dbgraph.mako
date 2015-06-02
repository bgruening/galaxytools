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

.unselectable {
    pointer-events: none;
}

.node {
    stroke: #fff;
    stroke-width: 1.5px;
}

.link {
    stroke: black;
}

.plink {
    stroke: black;

}

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
    text-align: right;
    padding: 3px;
    margin: 1px;
    color: white;
}


text {
    fill: white;
    font-weight: bolder;
    font-size: 10px;
    text-anchor: middle;
}

.background {
      fill: #fff;
}


svg {
    font: 10px sans-serif;
    background-color: rgb(200,200,200);
}

#dbgraph_sequence, #dbgraph_dotbracket {
	width:700px;font-family: courier;
}

</style>

<body>
<h1>Secondary Structure<h1>

<input type="text" title="sequence"   id="dbgraph_sequence"   readonly /><br />
<input type="text" title="dotbracket" id="dbgraph_dotbracket" readonly /><br />

<script type="text/javascript">

function os_to_dotbracket(os,sequence_length)
{
	var dotbracket = [];
	for(var i = 0; i < sequence_length; i++)
	{
		dotbracket.push('.');
	}
	
	for(var i = 0; i < os.length; i++)
	{
		if(os[i].source < os[i].target)
		{
			dotbracket[os[i]['source']] = '(';
			dotbracket[os[i]['target']] = ')';
		}
		else
		{
			dotbracket[os[i]['target']] = '(';
			dotbracket[os[i]['source']] = ')';
		}
	}
	
	return dotbracket.join('');
}

function parse_ps(img) {
    var maindic={},
        seqfound=false,
        sequence="",
        bpp=[],
        mlp=[];
        ps=[];
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
            a=f.split(" ")[0]-1;
            b=f.split(" ")[1]-1;
            var row = {};
            row[keys[0]]=parseInt(a);
            row[keys[1]]=parseInt(b);
            row[keys[2]]=1;
            mlp.push(row);
        }
    }
    console.log("source: " + mlp[0]["source"]);
    console.log("target: " + mlp[0]["target"]);
    for (var i=1; i < sequence.length; i++){
            var row = {};
            row["source"]=i-1;
            row["target"]=i;
            row["value"]=1;
            ps.push(row);
        }
    maindic["base-pairing-probabilities"]=bpp;
    maindic["optimal-structure"]=mlp;
    maindic["sequence"]=sequence;
    maindic["primary-structure"]=ps;
    
    return JSON.stringify(maindic);
}



$(document).ready(function() {
    var hdaId   = '${trans.security.encode_id( hda.id )}',
        hdaExt  = '${hda.ext}',
        ajaxUrl = "${h.url_for( controller='/datasets', action='index')}/" + hdaId + "/display?to_ext=" + hdaExt;


    %if hda.ext == 'rna_eps':
    $.ajax( ajaxUrl, {
        success: function(psImg)
        {
            var psjson=parse_ps(psImg);

        var bpm=JSON.parse(psjson);

        var os=bpm["optimal-structure"];
        var ps=bpm["primary-structure"];
       
        console.log(os.length);
        console.log(ps.length);

        var margin = {top: 40, right: 40, bottom: 10, left: 40},
            width = 700,
            height = 600;

        var force = d3.layout.force()
            .charge(-120)
            .linkDistance(12)
            .size([width, height]);
    
        var isequence=bpm.sequence.split("");
        var seq=isequence.map(function (d){return {"name":d}});
        
        $("#dbgraph_sequence").attr('value',bpm.sequence);
        $("#dbgraph_dotbracket").attr('value',os_to_dotbracket(os,bpm.sequence.length));
        
        var svg = d3.select("body").append("svg")
            .attr("width", width)
            .attr("height", height)
            .attr("style", "stroke-width:2")
            .attr("border", 3);
        var borderPath = svg.append("rect")
            .attr("x", 0)
            .attr("y", 0)
            .attr("height", height)
            .attr("width", width)
            .style("stroke", "black")
            .style("fill", "none")
            .style("stroke-width", "3");



        force
            .nodes(seq)
            .links(bpm["optimal-structure"].concat(bpm["primary-structure"]))
            .start();
            


        
        var link = svg.selectAll(".link")
            .data(bpm["optimal-structure"])
            .enter().append("line")
            .attr("stroke-dasharray","4,2")
            .attr("class", "link");

        var plink = svg.selectAll(".plink")
            .data(bpm["primary-structure"])
            .enter().append("line")
            .attr("class", "plink");

        var node = svg.selectAll(".node")
            .data(seq)
            .enter().append("circle")
            .attr("class", "node")
            .attr("r", 8)
            .style("fill", "red")
            .attr("cx",300)
            .attr("cy",300)
            .call(force.drag);

        var text = svg.selectAll(".text")
            .data(seq)
            .enter().append("text")
            .attr("class", "unselectable")
            .text(function(d){return ""+d.name;});

        node.append("title")
            .text(function(d) { return d.name; });


        force.on("tick", function() {
            link.attr("x1", function(d) { return d.source.x; })
                .attr("y1", function(d) { return d.source.y; })
                .attr("x2", function(d) { return d.target.x; })
                .attr("y2", function(d) { return d.target.y; });
            plink.attr("x1", function(d) { return d.source.x; })
                .attr("y1", function(d) { return d.source.y; })
                .attr("x2", function(d) { return d.target.x; })
                .attr("y2", function(d) { return d.target.y; });

            node.attr("cx", function(d) { return d.x; })
                .attr("cy", function(d) { return d.y; });
            text.attr("x", function(d) { return d.x; })
                .attr("y", function(d) { return d.y+4; });
        });
            
        }


	});
    %endif
});

</script>

</body>
</html>
