set -euo pipefail

# Start a MySQL 5.6 container with parameters matching the expectations of Cromwell's CI configuration:
# port, database name, user name, password etc.
docker run \
  -p 3306:3306 \
  --name cromwell_publish_db \
  -e MYSQL_ROOT_PASSWORD=test \
  -e MYSQL_DATABASE=cromwell_test \
  -e MYSQL_USER=cromwell \
  -e MYSQL_PASSWORD=test \
  -d mysql/mysql-server:5.6


