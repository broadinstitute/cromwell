#!/bin/bash

set -e

java $SYSTEMPROPS -jar $(find /cromwell | grep 'cromwell.*\.jar') server
