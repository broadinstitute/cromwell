## Persisting Data Between Restarts

### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Configuration Files](ConfigurationFiles.md)
* [Docker](https://docs.docker.com/engine/installation/)

### Goals

Cromwell remembers everything it knows!

### Let's get started!

- Start the MySQL docker container with the following line:

```bash
docker run -p 3306:3306 --name NameOfTheContainer -e MYSQL_ROOT_PASSWORD=YourPassword -e MYSQL_DATABASE=DatabaseName -e MYSQL_USER=ChooseAName -e MYSQL_PASSWORD=YourOtherPassword -d mysql/mysql-server:5.5
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

### Next steps

You might find the following tutorials interesting to tackle next:

* [Server Mode](ServerMode)
