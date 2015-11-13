#!/bin/bash

set -e
java $JAVA_OPTS -Djava.library.path=./native -javaagent:/cromwell/cromwell-*.jar -jar /cromwell/cromwell-*.jar server
