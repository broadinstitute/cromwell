CREATE USER 'cromwell'@'localhost' IDENTIFIED BY 'cromwell';
GRANT ALL PRIVILEGES ON cromwell_db.* TO 'cromwell'@'localhost' WITH GRANT OPTION;
CREATE USER 'cromwell'@'%' IDENTIFIED BY 'cromwell';
GRANT ALL PRIVILEGES ON cromwell_db.* TO 'cromwell'@'%' WITH GRANT OPTION;
CREATE USER 'cromwell'@'192.168.8.3' IDENTIFIED BY 'cromwell';
GRANT ALL PRIVILEGES ON cromwell_db.* TO 'cromwell'@'192.168.8.3' WITH GRANT OPTION;
