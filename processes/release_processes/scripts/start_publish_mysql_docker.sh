# Start a MySQL 5.6 container with parameters matching the expectations of Cromwell's CI configuration
# (port, database name, user name, password etc). Note this container is launched *without* `--detach`
# as that may lead to undesired behavior (the MySQL server could remain running in the background
# after a publish completes to produce an unwanted call cache hit for `publish_workflow.prepGithub`
# on a subsequent publish).
#
# This container should be terminated once the publish is complete and should clean up after itself as it
# is invoked with the `--rm` flag. This container appears to ignore `^C` but a combination of `^Z` and
# an appropriately parameterized `kill` does appear to be effective:
#
# 2021-04-14 12:22:57 0 [Note] mysqld (mysqld 5.6.51) starting as process 1 ...
# ^Z
# [1]+  Stopped                 ./start_publish_mysql_docker.sh
# scripts $ kill %1
# [1]+  Terminated: 15          ./start_publish_mysql_docker.sh
# scripts $

# Fail loudly if things go wrong https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425
set -euo pipefail

docker run \
  --publish 3306:3306 \
  --name cromwell_publish_db \
  --env MYSQL_ROOT_PASSWORD=test \
  --env MYSQL_DATABASE=cromwell_test \
  --env MYSQL_USER=cromwell \
  --env MYSQL_PASSWORD=test \
  --rm mysql/mysql-server:5.6
