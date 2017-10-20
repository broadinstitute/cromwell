#!/usr/bin/env bash

set -e
set -x

sbt clean assembly doc
