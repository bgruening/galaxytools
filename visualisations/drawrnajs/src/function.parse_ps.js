
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