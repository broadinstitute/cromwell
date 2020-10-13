
Cromwell can be thought of as performing three fundamental roles:

1. Front end - handling all REST requests
1. Runner - picking up and running workflows
1. Summarizer - summarizing metadata

In the simplest Cromwell deployment, a single Cromwell instance performs all three of these roles.
But it is also possible to create Cromwell deployments with many Cromwell instances, gaining advantages in
resiliency and scalability. There is only one restriction on roles: there must be exactly one Cromwell
instance performing the role of summarizer. Apart from this restriction, any instance in a multi-Cromwell
deployment can perform one or more of these roles. As a real world example, the current production Terra
deployment uses 7 Cromwell instances, each dedicated to a single role: 3 front ends, 3 runners, and one summarizer.

These roles are not explicit in Cromwell configuration, rather they are implied by configuration settings as
described below.

** Summarizer configuration **

The frequency of metadata summarization is determined by the value for the configuration key
`services.MetadataService.config.metadata-summary-refresh-interval`, which has a default value of `1 second`.
As stated above, there must be exactly one Cromwell instance performing the role of summarizer, so
all Cromwell instances which are not performing the summarizer role should specify `Inf` for this value.

** Runner configuration **

Cromwell instances in the runner role should periodically scan the workflow store to pick up and run
unclaimed workflows. The relevant configuration parameters are described below along with their default values.
The default values may be adequate for instances in the runner role, but will need to be overridden for
non-runner instances to effectively turn running off.

```hocon
system {
...
  # Number of seconds between polls of the workflow store.
  # Set this to a very large value for non-runners (e.g. 999999)
  new-workflow-poll-rate = 20

  # Cromwell will launch up to N submitted workflows at a time, regardless of how many open workflow slots exist
  # Set this to 0 for a non-runner.
  max-workflow-launch-count = 50

  # The maximum number of workflows to run concurrently.
  # Set this to 0 for a non-runner.
  max-concurrent-workflows = 5000
...
}
```

The documentation on [workflow heartbeats](https://cromwell.readthedocs.io/en/stable/Configuring/#workflow-heartbeats) describes how multiple Cromwell
runners collaborate to run workflows from a single workflow store.

** Front end configuration **

Cromwell instances should not require any configuration changes to operate in the front end role. If a particular
Cromwell instance is not intended to operate in the front end role then requests should not be directed to
that instance. If there are multiple front end instances then it may be desirable to configure a load balancer
in front of these instances to direct requests to the front end instances only.

