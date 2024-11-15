
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
                idx[name]['sequence'] += val;
            }
        }
        else if(state == 1)
        {
            if(val == "\n")
            {
                idx[name] = {'sequence':'','structure':''};
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
