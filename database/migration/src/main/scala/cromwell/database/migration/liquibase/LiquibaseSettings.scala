package cromwell.database.migration.liquibase

import liquibase.database.Database

case class LiquibaseSettings
(
  changeLogResourcePath: String,
  databaseChangeLogLockTableName: String = Database.databaseChangeLogLockTableName,
  databaseChangeLogTableName: String = Database.databaseChangeLogTableName
)
