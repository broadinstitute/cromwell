## Server

The default mode for most applications of Cromwell, suitable for production use. Server mode starts Cromwell as a web server that exposes REST endpoints. All features and APIs are available.

By default the server will be accessible at `http://localhost:8000`. See the [Server Section](Configuring#server) of the configuration for more information on how to configure it. A description of the endpoints can be found in the [API Section](api/RESTAPI).

Follow the [Server Tutorial](tutorials/ServerMode) to get your Cromwell server up and running in a few steps.

## Run

A good way to get started with Cromwell and experiment quickly. Run mode launches a single workflow from the command line and exits `0` or `1` to indicate the result. Appropriate for prototyping or demo use on a user's local machine. Features are limited and the web API is not supported.

Sending a `SIGINT` signal (via `CTRL-C` for example) will by default abort all running jobs and then exit.
This behavior can be configured, and is explained in more details in the [Abort](Configuring#abort) section of the configuration.
