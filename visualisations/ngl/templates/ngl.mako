<!DOCTYPE HTML>
<html>
<head>

<!--

NGL Viewer is a web application for molecular visualization. WebGL is employed to display molecules like proteins and 
DNA/RNA with a variety of representations. NGL is developed by Alexander Rose

To reference this tool please cite:

NGL Viewer: a web application for molecular visualization
10.1093/nar/gkv402

The source is MIT licensed and hosted under:
https://github.com/arose/ngl

Thanks Alex!

-->
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0">

    <title>${hda.name} | ${visualization_name}</title>
<%
    root = h.url_for( '/' )
%>

${h.stylesheet_link( root + 'plugins/visualizations/ngl/static/css/font-awesome.min.css' )}
${h.stylesheet_link( root + 'plugins/visualizations/ngl/static/css/main.css' )}
${h.stylesheet_link( root + 'plugins/visualizations/ngl/static/css/light.css', id='theme')}
${h.javascript_link( root + 'plugins/visualizations/ngl/static/js/ngl.full.min.js' )}

</head>

## ----------------------------------------------------------------------------
<body>
    <script>

        if ( ! Detector.webgl ) Detector.addGetWebGLMessage();

        document.addEventListener("DOMContentLoaded", function() {

            NGL.Stage.prototype.setTheme = function(){
                this.viewer.setBackground( "white" );
            }

            NGL.init( function(){

                var stage = new NGL.Stage();
                NGL.StageWidget( stage );

                stage.loadFile( 
                    window.location.protocol + '//' + window.location.host + 
                    "${h.url_for( controller='/datasets', action='index')}/${trans.security.encode_id( hda.id )}/display?to_ext=${hda.ext}",
                    undefined, undefined, undefined, 
                    {ext: "${hda.ext}", name: "${hda.name}"}
                );

            } );

        } );

    </script>
</body>
