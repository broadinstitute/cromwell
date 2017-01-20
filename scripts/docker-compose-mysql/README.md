# Cromwell server on MySQL Database

Uses docker-compose to link together a Cromwell docker image (built locally with `sbt docker` or available on [dockerhub](https://hub.docker.com/r/broadinstitute/cromwell/)) with a MySQL docker image.
To change the version of Cromwell used, change the tag in `compose/cromwell/Dockerfile`

## Local

`docker-compose up` from this directory will start a Cromwell server running on a mysql instance with local backend.
The default configuration file used can be found at `compose/cromwell/app-config/application.conf`
To override it, simply mount a volume containing your custom `application.conf` to `/app-config` (see `jes-cromwell/docker-compose.yml` for an example)

## JES

The `jes-cromwell` directory is an example of how to customize the original compose file with a configuration file and environment variables.
It uses the application default credentials of the host machine. To use it make sure your gcloud is up to date and that your [application-default credentials](https://developers.google.com/identity/protocols/application-default-credentials) are set up.
Then run `docker-compose -f docker-compose.yml -f jes-cromwell/docker-compose.yml up` to start a Cromwell server with a JES backend on MySQL.

## MySQL

The data directory in the MySQL container is mounted to `compose/mysql/data`, which allows the data to survive a `docker-compose down`.
To disable this feature, simply remove the `./compose/mysql/data:/var/lib/mysql` line in the volume section of `docker-compose.yml`.
Note that in such case, the data will still be preserved by a `docker-compose stop` that stops the container but doesn't delete it.

## Notes

To run Cromwell in the background, add `-d` at the end of the command:
`docker-compose up -d`

To then see the logs for a specific service, run `docker-compose logs -f <service>`. 
For example `docker-compose logs -f cromwell`.

For more information about docker compose: [Docker compose doc](https://docs.docker.com/compose/)
