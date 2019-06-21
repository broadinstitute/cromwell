#!/usr/bin/env bash

# Stores github credentials in a file
# https://git-scm.com/book/en/v2/Git-Tools-Credential-Storage

set -euo pipefail

program="$(basename "$0")"

usage () {
    cat <<USAGE
Usage: ${program} -t tokenFile -u user -e email [-n name] [-c credentialsFile]

Options:
    -t|--tokenFile <tokenFile>             : A file containing the GitHub token (required)
    -u|--user <user>                       : The GitHub user name (required)
    -e|--email <email>                     : The email used for commits (required)
    -n|--name <name>                       : The full name used for commits (optional)
    -c|--credentialsFile <credentialsFile> : The path to store the credentials (default ./githubCredentials)
USAGE
}

credentialsFile="$PWD/githubCredentials"


if ! OPTIONS=$(getopt -n "${program}" -o t:u:n:e:c: --long tokenFile:,user:,name:,email:,credentialsFile: -- "$@"); then
    usage
    exit 1
fi

eval set -- "${OPTIONS}"

while [[ $# -gt 0 ]]
do
    case "$1" in
        -t|--tokenFile       ) tokenFile="$2"; shift;;
        -u|--user            ) user="$2"; shift;;
        -e|--email           ) email="$2"; shift;;
        -n|--name            ) name="$2"; shift;;
        -c|--credentialsFile ) credentialsFile="$2"; shift;;
        --                   ) ;;
        *                    ) usage; exit 1;;
    esac
    shift
done

exit_error=1

if [[ -z ${tokenFile:+x} ]]; then
    echo "Error: Token file not specified" >&2
    exit_error=0
elif [[ ! -f "${tokenFile}" ]]; then
    echo "Error: Token file does not exist" >&2
    exit_error=0
fi

if [[ -z ${user:+x} ]]; then
    echo "Error: User not specified" >&2
    exit_error=0
fi

if [[ -z ${email:+x} ]]; then
    echo "Error: Email not specified" >&2
    exit_error=0
fi

if [[ ${exit_error} -eq 0 ]]; then
    echo
    usage
    exit 1
fi

echo "https://${user}:$(cat "${tokenFile}")@github.com" > "${credentialsFile}"
git config --global credential.helper "store --file ${credentialsFile}"
git config --global user.email "${email}"
if [[ -n ${name:+x} ]]; then
    git config --global user.name "${name}"
fi
