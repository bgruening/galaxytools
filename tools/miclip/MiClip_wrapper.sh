R --vanilla --slave --args $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} < /HOME/galaxy/galaxy-dist/tools/MiClip/MiClip.R > dump
if [ -f ${__tool_data_path__}clusters.csv ]
 then
   zip temp ${__tool_data_path__}log.txt
   zip temp ${__tool_data_path__}*.csv
   mv ${__tool_data_path__}temp.zip ${12}
 else 
   cat dump >&2
fi

