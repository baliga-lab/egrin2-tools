 #!/bin/bash          

 # initialize new composite database with
 # schema from individual cMonkey run
 # plus some additional tables and run_id
 # column

 # defaults to output in current working dir
 DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
 DBNAME="egrin2.db"
 echo "Initializing database $DBNAME @ $DIR"

 EGRIN2_DBS="$DIR/$(ls -d -R 2>/dev/null eco-ens*/)"
 #EGRIN_ROOT = $("eco-out-")
 #echo ${EGRIN2_DBS}