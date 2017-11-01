## Run

Run mode will run a single workflow from the command line, and exit when the workflow completes (successfully or not).
The exit code of the run command will be `0` if the workflow succeeds or is aborted, `1` if it fails.

Sending a `SIGINT` signal (via `CTRL-C` for example) will by default abort all running jobs and then exit.
This behavior can be configured, and is explained in more details in the [Abort](Configuring#abort) section of the configuration.

Run mode is a good way to get started with Cromwell and experiment quickly.
For more advanced use cases involving large workflows or multi-tenancy for example, Server mode is recommended. It also offers a larger set of features through its endpoints that are not available in run mode.

## Server

Server mode will start Cromwell as a web server that exposes REST endpoints.
By default the server will be accessible at `http://localhost:8000`. See the [Server Section](Configuring#server) of the configuration for more information on how to configure it.

A description of those endpoints can be found in the [API Section](api/RESTAPI).

Follow the [Server Tutorial](tutorials/ServerMode) to get your Cromwell server up and running in a few steps !
