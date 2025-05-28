## Persisting Data Between Restarts

### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Configuration Files](ConfigurationFiles.md)
* [Docker](https://docs.docker.com/engine/installation/)

### Let's get started!

- Start the MySQL docker container with the following line:

```bash
docker run -p 3306:3306 --name NameOfTheContainer -e MYSQL_ROOT_PASSWORD=YourPassword -e MYSQL_DATABASE=DatabaseName -e MYSQL_USER=ChooseAName -e MYSQL_PASSWORD=YourOtherPassword -d mysql:8.4
```

- Update your `application.conf` file.

```hocon
database {
  profile = "slick.jdbc.MySQLProfile$"
  db {
    driver = "com.mysql.cj.jdbc.Driver"
    url = "jdbc:mysql://localhost/DatabaseName?rewriteBatchedStatements=true&useSSL=false"
    user = "ChooseAName"
    password = "YourOtherPassword"
    connectionTimeout = 5000
  }
}
```
Add the line above, below the all other lines in your `application.conf`. Replace `"DatabaseName"`, `"ChooseAName"` and `"YourOtherPassword"` with the values you choose in step 2, preserving the double quotes.

Test it by running your server with the updated `application.conf`:
```bash
java -Dconfig.file=/path/to/application.conf/ -jar cromwell-[version].jar ...
```

### Rootless database option: Podman Quadlet

Podman uses [the same underlying technology](https://www.redhat.com/en/blog/containers-are-linux) as Docker, but its architecture allows running without root. Setting up Podman in your environment is beyond the scope of this tutorial, but once you have it going you can plug in this Quadlet.

[Learn about Quadlets in Podman Desktop](https://podman-desktop.io/blog/podman-quadlet).

Create this file in Podman Machine at path `/var/home/core/.config/containers/systemd/cromwell_database_3306.container`:

```
# cromwell_database_3306.container
[Container]
ContainerName=cromwell_database_3306
Environment=MYSQL_ROOT_PASSWORD=private MYSQL_USER=cromwell MYSQL_PASSWORD=test MYSQL_DATABASE=cromwell_test
Image=mirror.gcr.io/mysql:lts
PublishPort=3306:3306

# Use config from Cromwell repo
Volume=~/Projects/cromwell/src/ci/docker-compose/mysql-conf.d:/etc/mysql/conf.d

# Store data on the host for persistence and easy examination
Volume=~/Projects/local-cromwell-mysql:/var/lib/mysql

# Optional, periodically update image from the registry. Causes container restart.
# Requires one-time enablement on Podman Machine:
# $ systemctl enable podman-auto-update
# $ systemctl start podman-auto-update
AutoUpdate=registry

[Service]
Restart=always
```

### Next steps

You might find the following tutorials interesting to tackle next:

* [Server Mode](ServerMode)
