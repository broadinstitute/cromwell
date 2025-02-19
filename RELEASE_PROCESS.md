# Release Processes

## Core Process: How to Publish and Release Cromwell

This process is unrelated to Terra. Both GCP and Azure Terra update automatically. For more information see [here](https://support.terra.bio/hc/en-us/articles/9512163608731-Faster-Cromwell-updates-in-Terra-).

Manually cutting a numeric release is expected to continue being done on a cadence of every ~6 months. This is handled
through the Release Community Cromwell Github Action. Just run the action via the Github web UI - you must be a repo admin.

Before:
 * Communicate with the team in our private channel to ensure no one merges to `develop` while the release is ongoing.

After:
 * Announce the release in our public channel with context that it's for standalone Cromwells (the code is already in Terra).
 * The day after the release, confirm that [the Homebrew package](https://formulae.brew.sh/formula/cromwell) has the latest version. If it doesn't, start investigation by looking at [Homebrew PR's](https://github.com/Homebrew/homebrew-core/pulls?q=is%3Apr+cromwell).

## Bonus Processes

The Swagger client library is not part of our core publish/release process but can be performed from time to time, as required.

### How to Generate and Publish Swagger Client Library

**Note:** This part of publishing may or may not work for you until
[BT-38](https://broadworkbench.atlassian.net/browse/BT-38) is finished and this section is updated.

The first step is to generate the client library.  From the root of the repo run

```
./scripts/gen_java_client.sh
```

This generates the client library and runs the generated tests as well.  A successful run should end with something similar to

```
[debug] Test run finished: 0 failed, 0 ignored, 7 total, 0.007s
[info] Passed: Total 103, Failed 0, Errors 0, Passed 100, Skipped 3
[success] Total time: 4 s, completed Jun 26, 2019 3:01:19 PM
```

To publish to artifactory, first obtain the artifactory username and credentials. Credentials are located in Vault at path `secret/dsde/cromwell/common/cromwell-artifactory`. Then run

```
export ARTIFACTORY_USERNAME=<the-username>
export ARTIFACTORY_PASSWORD=<the-password>
./scripts/publish-client.sh
```

A SNAP version of the client library will be published in jFrog at [this path](https://broadinstitute.jfrog.io/ui/repos/tree/General/libs-release-local/org/broadinstitute/cromwell).
