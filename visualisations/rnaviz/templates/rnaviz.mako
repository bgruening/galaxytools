<!DOCTYPE html>
<html>
<head lang="en">
    <meta name="viewport" content="width=device-width,initial-scale=1,minimum-scale=1,maximum-scale=1"/>
    <title>Galaxy SQLite Data Viewer</title>

    <link rel="stylesheet" href="/plugins/visualizations/rnaviz/static/css/bootstrap.min.css">
    <link href="/plugins/visualizations/rnaviz/static/css/font-awesome.min.css" rel="stylesheet">
    <link rel="stylesheet" href="/plugins/visualizations/rnaviz/static/css/jquery.mCustomScrollbar.min.css">
    <link rel="stylesheet" href="/plugins/visualizations/rnaviz/static/css/rna.viz.css">
    
</head>
<body class="body-rna-viz">
    <div class="main-container">
        <div class="samples-overlay loader"><span> loading... </span></div>
    </div>
    <script src="/plugins/visualizations/rnaviz/static/js/lib/jquery-3.2.1.min.js"></script>
    <script src="/plugins/visualizations/rnaviz/static/js/lib/plotly.min.js"></script>
    <script src="/plugins/visualizations/rnaviz/static/js/lib/bootstrap.min.js"></script>
    <script src="/plugins/visualizations/rnaviz/static/js/lib/underscore.min.js"></script>
    <script src="/plugins/visualizations/rnaviz/static/js/lib/cytoscape.min.js"></script>
    <script src="/plugins/visualizations/rnaviz/static/js/lib/jquery.mCustomScrollbar.concat.min.js"></script>

    <script>
        $(document).ready(function () {
            var config = {
                href: document.location.origin,
                dataName: '${hda.name}',
                datasetID: '${trans.security.encode_id( hda.id )}',
                tableNames: {
                    % for table in hda.metadata.table_row_count:
                        "name": '${table}',
                    % endfor
                }
            };
            RNAInteractionViewer.loadData( config );
        });
    </script>
  
    <script src="/plugins/visualizations/rnaviz/static/js/visualize-alignment.js"></script>
    <script src="/plugins/visualizations/rnaviz/static/js/rna-viz.js"></script>

</body>
</html>
