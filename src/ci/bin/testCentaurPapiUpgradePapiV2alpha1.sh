#!/usr/bin/env bash

echo "$(tput setab 1)$(tput blink)BT-81: This sub-build was doing no useful work; it built Cromwell and launched Centaur but never found tests to run. However, tests for functionality implied by this sub-build's name may be desirable, see BT-81 for further information.$(tput sgr 0)"exit 0
