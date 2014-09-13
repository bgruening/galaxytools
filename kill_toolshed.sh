#!/bin/sh
rm -r /home/bag/projects/galaxy/galaxy-central/tool_deps/*
rm -r /home/bag/projects/galaxy/shed_tools/*
rm -rf /home/bag/projects/galaxy/galaxy-central/database/tmp/*

rm /home/bag/projects/galaxy/galaxy-central/integrated_tool_panel.xml

dbname="galaxy"
username="bag"
wherecond="bag"
psql $dbname $username << EOF
TRUNCATE tool_shed_repository CASCADE;
EOF


echo '<?xml version="1.0"?>
<toolbox tool_path="../shed_tools">
</toolbox>' > /home/bag/projects/galaxy/galaxy-central/shed_tool_conf.xml


