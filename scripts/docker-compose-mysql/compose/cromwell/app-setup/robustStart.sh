#! /bin/bash

CROMWELL_CMD=$1
RC=1

while [[ $RC -ne 0 ]]; do
  mysqladmin -h mysql-db ping > /dev/null 2>&1
  RC=$?
  if [[ $RC -eq 0 ]]; then
     echo "Connection established."
     break;
  else
    echo "MySQL server is not up yet. Waiting 5 seconds and retrying..."
    sleep 5
  fi
done

( 
  $CROMWELL_CMD
)
