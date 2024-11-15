<%
    # Use root for resource loading.
    root     = h.url_for( '/' )
    app_root = root + "plugins/visualizations/drawrnajs/static/"
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


${h.stylesheet_link( app_root + "drawrnajs/spectrum.css" )}
${h.stylesheet_link( app_root + "drawrnajs/style.css" )}

${h.javascript_link( app_root + 'drawrnajs/drawrnajs.0.3.5.min.js' )}

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


${h.javascript_link( app_root + 'function.os_to_dotbracket.min.js' )}
${h.javascript_link( app_root + 'function.parse_ps.min.js' )}
${h.javascript_link( app_root + 'function.redraw.min.js' )}
${h.javascript_link( app_root + 'function.draw_sequence_selector.min.js' )}
${h.javascript_link( app_root + 'function.get_fasta_index.min.js' )}
${h.javascript_link( app_root + 'function.get_dotbracket_index.min.js' )}
${h.javascript_link( app_root + 'function.get_connectivitytable_index.min.js' )}


<script type="text/javascript">

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
                
                var sequence = idx[Object.keys(idx)[0]]['sequence'];
                redraw(sequence,'.'.repeat(sequence.length));
            }
        });
    %elif hda.ext == 'dot-bracket' or hda.ext == 'dbn':
        $.ajax(ajaxUrl,
        {
            success: function(dotbracket_file)
            {
                var idx = get_dotbracket_index(dotbracket_file);
                draw_sequence_selector(idx);
                
                redraw(idx[Object.keys(idx)[0]]['sequence'],idx[Object.keys(idx)[0]]['structure']);
            }
        });
    %elif hda.ext == 'ct':
        $.ajax(ajaxUrl,
        {
            success: function(connectivitytable_file)
            {
                var idx = get_connectivitytable_index(connectivitytable_file);
                draw_sequence_selector(idx);
                
                redraw(idx[Object.keys(idx)[0]]['sequence'],idx[Object.keys(idx)[0]]['structure']);
            }
        });
    %endif
});

</script>

</body>
</html>
