#!/usr/bin/env bash

set -e
set -x

sbt +clean +package assembly +doc
