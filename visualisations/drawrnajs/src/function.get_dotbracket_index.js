
function get_dotbracket_index(dotbracket_content)
{
    // Gaurentees that the last sequence will be parsed even if no "\n" is the last char
    if(dotbracket_content[dotbracket_content.length-1] != "\n")
    {
        dotbracket_content = dotbracket_content + "\n";
    }
    
    idx = {};
    
    var state = 0;
    var name = '';
    
    for(var i = 0; i < dotbracket_content.length; i++)
    {
        var val = dotbracket_content[i];
        if(state == 0)
        {
            if(val == ">")
            {
                name = '';
            }
            else if(val != "\n")
            {
                name += val;
            }
            else
            {
                idx[name] = {'sequence':'','structure':''};
                state = 1
            }
        }
        else if(state == 1)
        {
            if(val != "\n")///@todo more sophisticated regex?
            {
                idx[name]['sequence'] += val;
            }
            else
            {
                state = 2;
            }
        }
        else if(state == 2)
        {
            if(val != "\n")
            {
                idx[name]['structure'] += val;
            }
            else
            {
                state = 0;
            }
        }
    }
    
    return idx;
}
