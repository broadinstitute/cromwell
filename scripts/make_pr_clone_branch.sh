#!/bin/bash

set -euo pipefail

export PR_TO_CLONE="$1"

git fetch -f origin "pull/${PR_TO_CLONE}/head:${PR_TO_CLONE}_pr_clone"
git checkout "${PR_TO_CLONE}_pr_clone"

echo "PR for hash $(rev-parse --verify HEAD)" > pr_salt
git add pr_salt
git commit -m "Salting PR branch"

echo "Created (or updated) pr-able branch for PR ${PR_TO_CLONE}."
echo "You can now force push this branch to run the CI suite for the changes in that PR."
echo "NOTE: Do NOT merge this branch into develop, merge the original"
echo

