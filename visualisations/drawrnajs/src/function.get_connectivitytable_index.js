
function get_connectivitytable_index(connectivitytable_content)
{
    // Guerentees that the last sequence will be parsed even if no "\n" is the last char
    if(connectivitytable_content[connectivitytable_content.length-1] != "\n")
    {
        connectivitytable_content = connectivitytable_content + "\n";
    }
    
    idx = {};
    
    var i = 0;
    var name = '';
    
    var lines = connectivitytable_content.match(/[^\r\n]+/g);
    
    for(var line_id in lines)
    {
        var line = lines[line_id];
        
        if(i == 0)// Header line
        {
            var columns = line.match(/([0-9]+)(?:\t|[ ]{2,})([^\n]+)/);
            name =  columns[2];
            i = parseInt(columns[1]);
            idx[name] = {'sequence':'','structure':'.'.repeat(i)};
        }
        else
        {
            var columns = line.match(/[^\t\ ]+/g);
            
            idx[name]['sequence'] += columns[1];
            
            var partner = parseInt(columns[4]);
            
            if(partner != 0)
            {
                var structure = idx[name]['structure'].split('');
                
                structure[partner - 1] = '(';
                structure[parseInt(columns[0]) - 1] = ')';
                
                idx[name]['structure'] = structure.join('');
            }
            
            i -= 1;
        }
    }
    
    return idx;
}
