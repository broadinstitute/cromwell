#!/bin/bash

set -e

java $JAVA_OPTS -jar $(find /cromwell | grep 'cromwell.*\.jar') server
