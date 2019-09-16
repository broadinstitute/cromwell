# How to Handle External Contributions

## Overview

- Decide whether the PR is adding community features or affects "supported" functionality.
- Review the concept
- Review the changes in the PR
- Run CI against the PR
- Cycle through Review and CI until satisfied
- Merge the PR

## Process

### Decide "Community" or "Supported"

**Community Supported:** If the PR only changes parts of Cromwell which are not part of the supported feature set, treat them with a little more
flexibility and with the assumption that the contributor probably knows better than us what they want the feature to do.

**Officially Supported:** If the PR changes core/supported features in Cromwell, review them as thoroughly as you would PRs from within the team. Remember
that one day you might need to support this code!

### Reviewing the Concept

Ask the questions:
 
- Will Cromwell be a better product with this change adopted. 
- Will it be better enough to warrant the time necessary to review the PR
  - Note: The answer to this is almost always a yes if the first answer was yes
  - However, overly long, opaque, or "risky" changes might benefit from requests to break the PR up and merge/review things in stages. 
  
### Review the changes in the PR

- For PRs changing "supported" features, treat it like any other PR coming from within the team.
  - Remember: we will have to support these changes in the future. Possibly forever!
- For PRs only making changes to features we don't officially support - be generous. But make sure:
  - That any new optional functionality is opt-in rather than becoming the default.
  - That any community features are flagged in documentation and config examples as originating from the community (and thus may not be supported by the team if bugs are found).
  
### Run CI against the PR

Note: inspired by the community answer [here](https://github.community/t5/How-to-use-Git-and-GitHub/Checkout-a-branch-from-a-fork/td-p/77).

- Problem: our CI will only run against branches of the `broadinstitute/cromwell` repo submitted by team members.
- To turn a community contribution into a PR that travis will run against:
  - Identify a reference to use for the remote branch and check it out. 
    - Example: for pull request 938 there will be a reference `refs/pull/938/head`
    - We can fetch that reference to a new branch using: `git fetch -f origin refs/pull/938/head:938_pr_clone`
  - Push the branch to github
  - Create a new PR for the clone branch. Indicate that it only exists to test a community contribution. 
    - Example title: `[PR 938 Clone] PR for CI only`

### Cycle through Review and CI

- If the community contributor makes changes following your reviews or the CI results:
  - Glance at the changes to make sure they still seem reasonable.
  - Make any additional comments
  - Re-fetch the remote reference: `git fetch -f origin refs/pull/860/head:938_pr_clone`
  - Push the changes back up to github to re-trigger the CI on your clone PR.
  
### Merge the PR

- Once the tests have completed successfully and the PR has two approvals, it can be merged.
- Remember to delete your branch clone PR (and the cloned branch itself too!)

## Shortcuts

### Git Command Shortcut

Note: also inspired by the community answer [here](https://github.community/t5/How-to-use-Git-and-GitHub/Checkout-a-branch-from-a-fork/td-p/77 and the reference gitconfig file [here](https://github.com/lee-dohm/dotfiles/blob/8d3c59004154571578c2b32df2cdebb013517630/gitconfig#L8)).

It's tedious to have to remember the syntax for `git fetch -f origin refs/pull/938/head:938_pr_clone` isn't it? Well
luckily you don't have to!

**Step 1:** add this line into your `~/.gitconfig` file under the `[alias]` section:
```
clone-pr = !sh -c 'git fetch -f origin pull/$1/head:$1_pr_clone && git checkout $1_pr_clone' -
```

**Step 2:** Your regular `git` command line now has new superpowers:
```
[develop] $ git clone-pr 938
From https://github.com/broadinstitute/cromwell
 * [new ref]             refs/pull/938/head -> 938_pr_clone
Switched to branch '938_pr_clone'

[938_pr_clone] $ git push
```

Note: The final `git push` command may not work as-is, depending on your ~/.gitconfig value of `push.default`. 
If it doesn't work then one of the following solutions may work:
  * Setting `git`'s `push.default` config value to be '`current`'.
  * Using `git push --set-upstream origin 5070_pr_clone` instead
  * Using `git push origin HEAD` instead

**Step 3:** If you need to re-sync your cloned PR against changes on their remote branch - no problem! The
exact same `git clone-pr 938` will *update* your local reference allowing you to push changes up to github easily!  
