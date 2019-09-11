#!/bin/sh
#   Use this script to ensure git-secrets are setup
#   https://minnie-kenny.readthedocs.io/

set -eu # -o pipefail isn't supported by POSIX

minnie_kenny_command_name=${0##*/}
minnie_kenny_quiet=0
minnie_kenny_strict=0
minnie_kenny_modify=0
minnie_kenny_gitconfig="minnie-kenny.gitconfig"

usage() {
  if [ ${minnie_kenny_quiet} -ne 1 ]; then
    cat <<USAGE >&2
Usage:
    ${minnie_kenny_command_name}
    -f | --force                Modify the git config to run git secrets
    -n | --no-force             Do not modify the git config, only verify installation
    -s | --strict               Require git-secrets to be setup or fail
    -q | --quiet                Do not output any status messages
    -i | --include=FILE         Path to the include for git-config (default: "minnie-kenny.gitconfig")
USAGE
  fi
  exit 1
}

run_command() { if [ ${minnie_kenny_quiet} -ne 1 ]; then "$@"; else "$@" >/dev/null 2>&1; fi; }
echo_out() { if [ ${minnie_kenny_quiet} -ne 1 ]; then echo "$@"; fi; }
echo_err() { if [ ${minnie_kenny_quiet} -ne 1 ]; then echo "$@" 1>&2; fi; }

process_arguments() {
  while [ $# -gt 0 ]; do
    case "$1" in
      -q | --quiet)
        minnie_kenny_quiet=1
        shift 1
        ;;
      -s | --strict)
        minnie_kenny_strict=1
        shift 1
        ;;
      -f | --force)
        minnie_kenny_modify=1
        shift 1
        ;;
      -n | --no-force)
        minnie_kenny_modify=0
        shift 1
        ;;
      -i)
        shift 1
        minnie_kenny_gitconfig="${1:-}"
        if [ "${minnie_kenny_gitconfig}" = "" ]; then break; fi
        shift 1
        ;;
      --include=*)
        minnie_kenny_gitconfig="${1#*=}"
        shift 1
        ;;
      --help)
        usage
        ;;
      *)
        echo_err "Unknown argument: $1"
        usage
        ;;
    esac
  done

  if [ "${minnie_kenny_gitconfig}" = "" ]; then
    echo_err "Error: you need to provide a git-config include file."
    usage
  fi
}

# Exits if this system or directory is not setup to run git / minnie-kenny.sh / git-secrets
validate_setup() {
  if ! command -v git >/dev/null 2>&1; then
    if [ ${minnie_kenny_strict} -eq 0 ]; then
      echo_out "\`git\` not found. Not checking for git-secrets."
      exit 0
    else
      echo_err "Error: \`git\` not found."
      exit 1
    fi
  fi

  minnie_kenny_is_work_tree="$(git rev-parse --is-inside-work-tree 2>/dev/null || echo false)"

  if [ "${minnie_kenny_is_work_tree}" != "true" ]; then
    if [ ${minnie_kenny_strict} -eq 0 ]; then
      echo_out "Not a git working tree. Not checking for git-secrets."
      exit 0
    else
      echo_err "Error: Not a git working tree."
      exit 1
    fi
  fi

  # Get the git absolute directory, even when older versions of git do not support --absolute-git-dir
  minnie_kenny_git_dir="$(cd "$(git rev-parse --git-dir)" && pwd)"
  if [ ! -f "${minnie_kenny_git_dir}/../${minnie_kenny_gitconfig}" ]; then
    echo_err "Error: ${minnie_kenny_gitconfig} was not found next to the directory ${minnie_kenny_git_dir}"
    exit 1
  fi

  if ! command -v git-secrets >/dev/null 2>&1; then
    echo_err "\`git-secrets\` was not found while \`git\` was found." \
      "\`git-secrets\` must be installed first before using ${minnie_kenny_command_name}." \
      "See https://github.com/awslabs/git-secrets#installing-git-secrets"
    exit 1
  fi
}

# Echo 1 if the hook contains a line that starts with "git secrets" otherwise echo 0
check_hook() {
  path="${minnie_kenny_git_dir}/hooks/$1"
  if grep -q "^git secrets " "${path}" 2>/dev/null; then
    echo 1
  else
    echo 0
  fi
}

# Ensures git secrets hooks are installed along with the configuration to read in the minnie-kenny.gitconfig
check_installation() {
  expected_hooks=0
  actual_hooks=0

  for path in "commit-msg" "pre-commit" "prepare-commit-msg"; do
    increment=$(check_hook ${path})
    actual_hooks=$((actual_hooks + increment))
    expected_hooks=$((expected_hooks + 1))
  done

  if [ 0 -lt ${actual_hooks} ] && [ ${actual_hooks} -lt ${expected_hooks} ]; then
    # Only some of the hooks are setup, meaning someone updated the hook files in an unexpected way.
    # Warn and exit as we cannot fix this with a simple `git secrets --install`.
    echo_err "Error: git-secrets is not installed into all of the expected git hooks." \
      "Double check the 'commit-msg' 'pre-commit' and 'prepare-commit-msg' under the directory" \
      "${minnie_kenny_git_dir}/hooks and consider running \`git secrets --install --force\`."
    exit 1
  fi

  # Begin checking for fixable errors
  found_fixable_errors=0

  if [ ${actual_hooks} -eq 0 ]; then
    if [ ${minnie_kenny_modify} -eq 1 ]; then
      run_command git secrets --install
    else
      echo_err "Error: git-secrets is not installed into the expected git hooks" \
        "'commit-msg' 'pre-commit' and 'prepare-commit-msg'."
      found_fixable_errors=1
    fi
  fi

  # Allow the minnie-kenny.gitconfig in `git secrets --scan`
  if ! git config --get-all secrets.allowed | grep -Fxq "^${minnie_kenny_gitconfig}:[0-9]+:"; then
    if [ ${minnie_kenny_modify} -eq 1 ]; then
      run_command git config --add secrets.allowed "^${minnie_kenny_gitconfig}:[0-9]+:"
    else
      echo_err "Error: The expression '^${minnie_kenny_gitconfig}:[0-9]+:' should be allowed by git secrets."
      found_fixable_errors=1
    fi
  fi

  # Allow minnie-kenny.gitconfig to appear in `git secrets --scan-history`
  if ! git config --get-all secrets.allowed | grep -Fxq "^[0-9a-f]+:${minnie_kenny_gitconfig}:[0-9]+:"; then
    if [ ${minnie_kenny_modify} -eq 1 ]; then
      run_command git config --add secrets.allowed "^[0-9a-f]+:${minnie_kenny_gitconfig}:[0-9]+:"
    else
      echo_err "Error: The expression '^[0-9a-f]+:${minnie_kenny_gitconfig}:[0-9]+:' should be allowed by git secrets."
      found_fixable_errors=1
    fi
  fi

  if ! git config --get-all include.path | grep -Fxq "../${minnie_kenny_gitconfig}"; then
    if [ ${minnie_kenny_modify} -eq 1 ]; then
      run_command git config --add include.path "../${minnie_kenny_gitconfig}"
    else
      echo_err "Error: The path '../${minnie_kenny_gitconfig}' should be an included path in the git config."
      found_fixable_errors=1
    fi
  fi

  if [ ${found_fixable_errors} -ne 0 ]; then
    echo_err "Error: The above errors may be fixed by re-running ${minnie_kenny_command_name} with -f / --force."
    exit 1
  fi
}

main() {
  process_arguments "$@"
  validate_setup
  check_installation
}

main "$@"
