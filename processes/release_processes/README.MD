# Release Processes

## How to update these processes

Have a better idea about how the deployment processes should work? 
See our "updating the process" [process](../README.MD)!

## Core Process: How to Publish and Release Cromwell

If you're the lucky "Release Rotation Engineer" this time, you should do the following three things:

First: Run the publish script to create a new version of Cromwell.
Second: Run through the "How to Release Cromwell into Firecloud" process.
Third: Run through the "How to Deploy Cromwell in CAAS prod" process.

### How to publish a new Cromwell version

#### Get a Github token

The release WDL uses a github token to perform actions on your behalf.
* Follow the directions [here](https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line).
    * When asked, give the token the full set of `repo` permissions and also one of either `user` or `user:email`.
* Copy the token and save it to a file on your local machine. You'll use it in a second.

#### Prepare a temporary `release` directory

This is optional, but I find it useful. Make or copy the following files into some temporary `releases/` directory:

* A cromwell jar file, preferably the most recent Cromwell version.
* A copy of the workflow file to run (https://github.com/broadinstitute/cromwell/blob/develop/publish/publish_workflow.wdl)
* An inputs json like this:

```json
{
  "publish_workflow.githubToken": "<<GITHUB TOKEN VALUE>>",
  "publish_workflow.majorRelease": true,
  "publish_workflow.publishHomebrew": true,
  "publish_workflow.publishDocker": "broadinstitute/cromwell-publish:latest",
  "publish_workflow.organization": "broadinstitute"
}
```

#### Make sure docker will have enough memory

I had to follow the instructions [here](https://docs.docker.com/docker-for-mac/#resources) to increase my docker memory.
I chose to increase it from 2GB to 8GB but 4 might have been sufficient.

#### Run the `publish_workflow.wdl` Workflow

* Run the Cromwell instance using the local-backend.
    * In server mode with a persistent backing database is probably a good idea - it will allow call caching to happen if you need to restart for any reason.
* Submit the workflow to Cromwell along with the inputs file.

#### Make sure it all went swimmingly

* Check that the workflow succeeded
* Check that there's now a new Cromwell release listed [here](https://github.com/broadinstitute/cromwell/releases).
* Check that there's a homebrew PR for the new Cromwell version [here](https://github.com/Homebrew/homebrew-core/pulls) (and that it passes their CI)  
* Have a look at https://travis-ci.com/broadinstitute/cromwell/branches. There should be a build for the new release tag that publishes new docker images for the version.

### How to Release Cromwell into Firecloud

**Note:** If the Cromwell CHANGELOG indicates that the upgrade might take some time (eg because of a database migration), checking in with the release engineer
and user support/comms to let them know that the upgrade may involve downtime is also required. You may need to help draft an impact statement and co-ordinate timing
the deploy to make sure user impact in minimized.

**Note:** How to accomplish some of these steps might be non-obvious to you. If so, refer to the additional details in the [full document](https://docs.google.com/document/d/1EEzwemE8IedCplIwL506fiqXr0262Pz4G0x6Cr6V-5E). 

![firecloud-develop](firecloud-develop.dot.png) 

#### Performance Testing

This performance testing process happens as part of the release into FireCloud.

* **Note:** This process is non-blocking and optional for HOTFIX releases.

1. Look at the [performance testing spreadsheet](https://docs.google.com/spreadsheets/d/1ksSHJdODuqqmvJ_94zZdc33Ai1Xv_53XzpsDd8jdmB0/edit#gid=479925276)
for any performance changes since the previous release.

2. A small amount of jitter is acceptable. So long as the lines are holding stable (or improving!), we consider Cromwell performance to not have regressed.

    * **Note:** If the spreadsheet is out of date, go to the testResults.tsv artifacts for recent runs of [Cromwell Perf: Limit Finder](https://fc-jenkins.dsp-techops.broadinstitute.org/job/cromwell-perf-limit-finder/)
runs and fill in any missing entries in the table.
 

### How to Deploy Cromwell in CAAS prod

**Note:** If the Cromwell CHANGELOG indicates that the upgrade might take some time (eg because of a database migration), checking in with the CaaS users
to let them know that the upgrade is about to happen is a good idea.

Deploying to CAAS is detailed in the [Quick CAAS Deployment Guide](https://docs.google.com/document/d/1s0YC-oohJ7o-OGcgnH_-YBtIEKmLIPTRpG36yvWxUpE) 

## Bonus Processes

The swagger client library is not part of our core publish/release process but can be performed from time to time, as required.

### How to Generate and Publish Swagger Client Library

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

To publish to artifactory, first obtain the artifactory username and credentials.  Then run

```
export ARTIFACTORY_USERNAME=<the-username>
export ARTIFACTORY_PASSWORD=<the-password>
./scripts/publish-client.sh
```
