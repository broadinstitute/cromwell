## Persisting Data Between Restarts

### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Configuration Files](ConfigurationFiles.md)

### Goals

At the end of this tutorial you'll have added a database to your configuration file so that Cromwell's knoweldge of workflow state and metadata for previous workflow executions is persisted even between invocations of Cromwell.

### Let's get started!


### Community-contributed 

A short guide to run MySQL-DB in a docker-container and use it for cromwell-server. 

Prerequisite:

* Cromwell >= 29 (https://github.com/broadinstitute/cromwell/releases)
* docker already installed (https://docs.docker.com/engine/installation/)

1. Pull the MySQL docker image from dockerhub:
`docker pull mysql:5.5`

2. Start the MySQL docker container with the following line:
```
$ docker run -p 3306:3306 --name NameOfTheContainer -e MYSQL_ROOT_PASSWORD=YourPassword-e MYSQL_DATABASE=DatabaseName -e MYSQL_USER=ChooseAName -e MYPASSWORD=YourOtherPassword -d mysql/mysql-server:latest
```
The `-p` option maps the standard port 3306 from the container to the 3306 on your machine running the container, this is necessary for Cromwell to communicate with the database inside the docker container. You can choose names and passwords freely but keep them in mind for the cromwell configuration file.

3. Update your `application.conf` file.

```
database { profile = "slick.jdbc.MySQLProfile$" db{ driver = "com.mysql.jdbc.Driver" url = "jdbc:mysql://localhost/DatabaseName?useSSL=false" user = "ChooseAName" password = "YourOtherPassword" connectionTimeout = 5000 } }
```
Add the line above, below the all other lines in your application.conf. Replace "DatabaseName","ChooseAName" and "YourOtherPassword" with the values you choose in step 2, preserving the double qoutes.

Test it, by running your server with the updated application.conf.
```
$ java -Dconfig.file=/path/to/application.conf/ -jar ...
```

The guide was tested on Debian 8.


### Easy trick to make DB persistent


### Next steps

You might find the following tutorials interesting to tackle next:

* [Server Mode](ServerMode.md)
