# Release Processes

## How to update these processes

Have a better idea about how the deployment processes should work? 
See our "updating the process" [process](../README.MD)!

## Core Process: How to Publish and Release Cromwell

If you're the lucky "Release Rotation Engineer" this time, you should do the following four things:

1. Check that all most recent executions of Cromwell's Jenkins tests have succeeded, and that all those executions were 
against the most recent commit in "develop" branch. There are 2 groups of Jenkins tests which have to be checked (you 
will need to be on the Broad internal network or VPN to open the following links):
    1. Tests for various backends supported by Cromwell. Log into Jenkins [here](https://fc-jenkins.dsp-techops.broadinstitute.org), check the tests [here](https://fc-jenkins.dsp-techops.broadinstitute.org/job/cromwell-cron-parent/).
    1. Tests for Cromwell in Terra environment. Log into Jenkins [here](https://fc-jenkins.dsp-techops.broadinstitute.org), check the tests [here](https://fc-jenkins.dsp-techops.broadinstitute.org/view/Batch/).
1. [Run the publish script to create a new version of Cromwell](#how-to-publish-a-new-cromwell-version)  
1. [Run through the "How to Release Cromwell into Firecloud/Terra" process](#how-to-release-cromwell-into-firecloud--terra)  
1. [Run through the "How to Deploy Cromwell in CaaS prod" process](#how-to-deploy-cromwell-in-caas-staging-and-caas-prod)  

### How to publish a new Cromwell version

#### Announce your intentions

Post a message in `#ftfy-private` letting people know that a publish is imminent in case there are PRs they want to get
wrapped up and merged to `develop` to go out in the forthcoming version of Cromwell.

#### Get a Github token

The release WDL uses a github token to perform actions on your behalf.
* Follow the directions [here](https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line).
  When asked give the token these exact permissions ([example screenshot](github_token_scopes.png)):
  * the full set of `repo` permissions
  * the `workflow` permission
  * only the `read:org` permission (do not let this token administer your organizations!)
  * only the `user:email` permission (do not let this token administer your organizations!)
* Copy the token and save it to a file on your local machine. You'll use it in a second.

#### Prepare a temporary `release` directory

This is optional, but I find it useful. Make or copy the following files into some temporary `releases/` directory:

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

#### Make sure Docker will have enough memory

Follow the instructions [here](https://docs.docker.com/docker-for-mac/#resources) to increase Docker memory.
Ensure you have at least 8GB; 4GB is not sufficient.

#### Let people know the publish is underway

Post another message in `#ftfy-private` that the release is underway, asking everyone to hold off merges to develop until
the release is published.

#### Run the `publish_workflow.wdl` Workflow

Run Cromwell in server mode with a persistent backing database, using Docker containers. This allows call caching to happen if you need to restart for any reason. 
See instructions for using a Dockerized MySQL server and CI config [here](#cromwell-setup-for-publishing).

Note that the publish workflow is quite resource intensive; it's a good idea to
shut down other resource intensive apps before launching it to avoid painfully slow or failed executions.

Using the Swagger API, submit the workflow to Cromwell along with the inputs file.

#### Make sure it all went swimmingly

* Check that the workflow succeeded
* Check that there's now a new Cromwell release listed [here](https://github.com/broadinstitute/cromwell/releases).
* If publishing to Homebrew, check that there's a Homebrew PR for the new Cromwell version [here](https://github.com/Homebrew/homebrew-core/pulls) (and that it passes their CI)
* Look [in Travis](https://app.travis-ci.com/github/broadinstitute/cromwell/branches) for the release tag build that will publish Docker images for the new version.
* Let `#ftfy-private` know that it's okay to resume merges to `develop`.
* Announce release in `#dsp-batch`, set expectations about when the new version will be available in Terra.

### How to Release Cromwell into Firecloud / Terra

**Note:** If the Cromwell CHANGELOG indicates that the upgrade might take some time (e.g., because of a database migration), checking in with the release engineer
and user support/comms to let them know that the upgrade may involve downtime is also required. You may need to help draft an impact statement and co-ordinate timing
the deploy to make sure user impact in minimized.

**Note:** How to accomplish some of these steps might be non-obvious to you (e.g., generating the release notes).
If so, refer to the additional details in the [full document](https://docs.google.com/document/d/1EEzwemE8IedCplIwL506fiqXr0262Pz4G0x6Cr6V-5E). 

![firecloud-develop](firecloud-develop.dot.png) 

### How to Deploy Cromwell in CaaS staging and CaaS prod

CaaS is "Cromwell as a Service". It is used by a couple of Broad teams (Pipelines and Epigenomics), though the long-term plan is for those teams to migrate to using Terra.

**Note:** If the Cromwell CHANGELOG indicates that the upgrade might take some time (e.g., because of a database migration), checking in with the CaaS users
to let them know that the upgrade is about to happen is a good idea.

Deploying to CaaS is detailed in the [Quick CaaS Deployment Guide](https://docs.google.com/document/d/1s0YC-oohJ7o-OGcgnH_-YBtIEKmLIPTRpG36yvWxUpE) 

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

To publish to artifactory, first obtain the artifactory username and credentials.  Then run

```
export ARTIFACTORY_USERNAME=<the-username>
export ARTIFACTORY_PASSWORD=<the-password>
./scripts/publish-client.sh
```

## Cromwell setup for publishing

One can run a publishing-friendly Cromwell using a containerized MySQL server and config files
from Cromwell's CI. See the `start_publish_mysql_docker.sh` and `start_publish_cromwell.sh` 
scripts under `release_processes/scripts`.
