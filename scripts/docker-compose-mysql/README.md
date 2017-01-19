# Cromwell server on MySQL Database

Uses docker-compose to link together a cromwell docker image (built locally with `sbt docker` or available on [dockerhub](https://hub.docker.com/r/broadinstitute/cromwell/)) with a MySQL docker image.
To change the version of cromwell used, change the tag in `compose/cromwell/Dockerfile`

## Local

`docker-compose up` from this directory will start a cromwell server running on a mysql instance with local backend.
The default configuration file used can be found at `compose/cromwell/app-config/application.conf`
To override it, simply mount a volume containing your custom `application.conf` to `/app-config` (see the Jes section below)

## Jes

The `jes-cromwell` directory is an example of how to customize the original compose file with a configuration file and environment variables.
It uses the application default credentials of the host machine. To use it make sure your gcloud is up to date and that your [application-default credentials](https://developers.google.com/identity/protocols/application-default-credentials) are setup.
Then run `docker-compose -f docker-compose.yml -f jes-cromwell/docker-compose.yml up` to start a cromwell server with a Jes backend on MySQL.

## Notes

To run cromwell in the background, add `-d` at the end of the command:
`docker-compose up -d`

To then see the logs for a specific service, run `docker-compose logs <service> -f`. 
For example `docker-compose logs cromwell -f`.

For more information about docker compose: [Docker compose doc](https://docs.docker.com/compose/)
