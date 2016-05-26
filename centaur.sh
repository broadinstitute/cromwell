#!/bin/bash

CMD='sbt'
if [ -n "$1" ]; then
    CMD=$1
fi

docker run -it --rm \
    -e "JAVA_TOOL_OPTIONS=-Dfile.encoding=UTF8i -Dcentaur.cromwellUrl=http://cromwell:8000" \
    --link centaur_cromwell_1:cromwell \
    broadinstitute/centaur:dev \
    $CMD
