
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
    runButton.addEventListener('click', function()
    { 
        document.getElementById('ALERT').value = "";
        var input = rna.io.getInputSequences();
        if(rna.io.checkConditions(input)){
            struct = rna.t.transformDotBracket(input[0], input[1]);
            app({graph: struct, el: cy, doc: document, win: window});
        }
    }, false);
}