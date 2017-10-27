Cromwell provides two ways to run workflows.

## Run

Run mode will run a single workflow, and exit when the workflow completes (successfully or not).
The exit code of the run command will be `0` if the workflow succeeds or is aborted, `1` if it fails.

Sending a `SIGINT` signal (via `CTRL-C` for example) will by default abort all running jobs and then exit.
This behavior can be configured, and is explained in more details in the [Abort](Configuring#Abort) section of the configuration.