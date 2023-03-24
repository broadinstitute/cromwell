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

- Opening a PR triggers CI to run automatically. However, if the PR was not opened by a trusted contributor, it will wait in a Github Actions queue until a trusted contributor explicitly allows CI tests run.

### Cycle through Review and CI

- If the community contributor makes changes following your reviews or the CI results:
  - Glance at the changes to make sure they still seem reasonable.
  - Make any additional comments
  - Rerun CI.
  
### Merge the PR

- Once the tests have completed successfully and the PR has two approvals, it can be merged.
- Remember to delete your branch clone PR (and the cloned branch itself too!)

