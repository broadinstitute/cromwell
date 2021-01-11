#!/usr/bin/env bash

echo "$(tput setab 1)$(tput blink)BT-84: This sub-build currently does no useful work; it builds Cromwell and launches Centaur but finds no matching tests to run. See the linked ticket for further information.$(tput sgr 0)"
exit 0
