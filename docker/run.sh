#!/bin/bash

set -e

java -jar $(find /cromwell | grep 'cromwell.*\.jar') server
