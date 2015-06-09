<%
    # Use root for resource loading.
    root                      = h.url_for( '/' )
    app_root                  = root + "plugins/visualizations/drawrnajs/static/"
%>
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


<meta charset="UTF-8">
<link type="text/css" rel="stylesheet" href="http://workmen.biojs.net/github/bene200/drawrnajs/master/style/style.css">
<link type="text/css" rel="stylesheet" href="http://workmen.biojs.net/github/bene200/drawrnajs/master/style/spectrum.css">

<script src="https://wzrd.in/bundle/drawrnajs@0.3.4"></script>


<!--
example page:
http://workmen.biojs.net/demo/drawrnajs/simple
-->

</head>
<body>
	<div class="input">
		<input type="text" value="" class="alertfield" id="ALERT" readonly>
		<input class="textbox" id="SEQ_BOX">
		<input class="textbox" id="DOTBR_BOX">  
	</div>
	<div class="main">
		<div class="container">
			<h2>Options</h1>
			<ul> 
				<li> A  <input type="text" id="acolor"> </li>
				<li> C  <input type="text" id="ccolor"> </li>
				<li> G  <input type="text" id="gcolor"> </li>
				<li> U  <input type="text" id="ucolor"> </li>
			</ul>
			<p>To add a hydrogen bond, first press "N". Then click on the nucleotides from which to which you want to create the bond. </p>
			<p>The Dot-Bracket Notation will be changed to fit the new structure automatically.</p>
			
			<input class="button" id="PERFORM_VIS" name="Perform Vis" value="Show structure">
			<input class="button" id="CENTER" name="CenterView" value="Center view">
			<input class="button" id="EXPORT" name="ExportStructure" value="Export as PNG">

			<div>
				<br>
				Color of selected nucleic acids <input type="text" id="selcolor"><br><br>
				<input class="button" id="SELCOL" name="Change" value="Change">
			</div>
		</div>
		<div class="drawfield">
			<div id="cy"></div>
		</div>
	</div>



<script type="text/javascript">
///@todo make it static in a separate .js file and run grunt over it
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

///@note this function is taken from:
///@url https://github.com/bgruening/galaxytools/tree/master/visualisations/dbgraph
function parse_ps(img)
{
	var maindic={},
		seqfound=false,
		sequence="",
		bpp=[],
		mlp=[];
		ps=[];
	var keys=["source","target","value"];
	
	var dic= img.split("\n");
	for (var i=0;i<dic.length ; i++)
	{
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

function redraw(sequence,structure)
{
	var rna = require("drawrnajs");
	var app = rna.vis;
	
	while(structure.length < sequence.length)
	{
		structure += '.';
	}
	
	document.getElementById('SEQ_BOX').value = sequence;
	document.getElementById('DOTBR_BOX').value = structure;
	
	//init colors
	$("#acolor").spectrum({ color: "#64F73F" });
	$("#ccolor").spectrum({ color: "#FFB340" });
	$("#gcolor").spectrum({ color: "#EB413C" });
	$("#ucolor").spectrum({ color: "#3C88EE" });
	$("#selcolor").spectrum({color: "#F6F6F6"});

	//init alert box
	document.getElementById('ALERT').value = "";

	var input = rna.io.getInputSequences();
	var struct = rna.t.transformDotBracket(input[0], input[1]);

	var cy = document.getElementById('cy');
	cy.style.width = "60%";

	app({graph: struct, el: cy, doc: document, win: window});

	var runButton = document.getElementById('PERFORM_VIS');
	runButton.readOnly = true;
	runButton.addEventListener('click', function(){ 
		document.getElementById('ALERT').value = "";
		var input = rna.io.getInputSequences();
		if(rna.io.checkConditions(input)){
			struct = rna.t.transformDotBracket(input[0], input[1]);
			app({graph: struct, el: cy, doc: document, win: window});
		}
	}, false);
}

function draw_sequence_selector(sequence_index)
{
	var keys = Object.keys(sequence_index);
	var n = keys.length;
	
	///@todo use $().data()
	var select = $('<select class="textbox" onchange="redraw(idx[Object.keys(idx)[$(this).val()]],\'\');" />');
	
	for(var i = 0; i < 2; i++)
	{
		var option = $('<option />');
		option.attr('value',i);
		option.text((i+1)+'. '+keys[i]);
		
		select.append(option);
	}
	
	$( ".input" ).prepend(select);
}

function get_fasta_index(fasta_content)
{
	// Gaurentees that the last sequence will be parsed even if no "\n" is the last char
	if(fasta_content[fasta_content.length-1] != "\n")
	{
		fasta_content = fasta_content + "\n";
	}
	
	idx = {};
	
	var state = 2;
	var name = '';
	
	for(var i = 0; i < fasta_content.length; i++)
	{
		var val = fasta_content[i];
		if(state == 2)
		{
			if(val == ">")
			{
				name = '';
				state = 1;
			}
			else if(val != "\t" && val != ' ' && val != "\n")
			{
				idx[name] += val;
			}
		}
		else if(state == 1)
		{
			if(val == "\n")
			{
				idx[name] = '';
				state = 2;
			}
			else
			{
				name += val;
			}
		}
	}
	
	return idx;
}

$(document).ready(function()
{
	var hdaId   = '${trans.security.encode_id( hda.id )}',
		hdaExt  = '${hda.ext}',
		ajaxUrl = "${h.url_for( controller='/datasets', action='index')}/" + hdaId + "/display?to_ext=" + hdaExt;
	
	%if hda.ext == 'rna_eps':
		$.ajax(ajaxUrl,
		{
			success: function(psImg)
			{
				var psjson=parse_ps(psImg);
				var bpm=JSON.parse(psjson);

				var os=bpm["optimal-structure"];
				var ps=bpm["primary-structure"];
				
				redraw(bpm.sequence,os_to_dotbracket(os,bpm.sequence.length));
			}
		});
	%elif hda.ext == 'fasta':
		$.ajax(ajaxUrl,
		{
			success: function(fasta_file)
			{
				var idx = get_fasta_index(fasta_file);
				draw_sequence_selector(idx);
				
				redraw(idx[Object.keys(idx)[0]],'');
			}
		});
	%endif
});

</script>

</body>
</html>
