
function draw_sequence_selector(sequence_index)
{
    var keys = Object.keys(sequence_index);
    var n = keys.length;
    
    ///@todo use $().data()
    var select = $('<select class="textbox" onchange="redraw(idx[Object.keys(idx)[$(this).val()]][\'sequence\'],idx[Object.keys(idx)[$(this).val()]][\'structure\']);" />');
    
    for(var i = 0; i < n; i++)
    {
        var option = $('<option />');
        option.attr('value',i);
        option.text((i+1)+'. '+keys[i]);
        
        select.append(option);
    }
    
    $( ".input" ).prepend(select);
}