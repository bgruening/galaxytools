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