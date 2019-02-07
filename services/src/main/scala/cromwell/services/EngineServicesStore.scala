package cromwell.services

import cromwell.database.migration.liquibase.LiquibaseSettings
import cromwell.database.slick.EngineSlickDatabase
import cromwell.database.sql.EngineSqlDatabase

object EngineServicesStore {
  lazy val EngineLiquibaseSettings = LiquibaseSettings("changelog.xml")

  import ServicesStore.EnhancedSqlDatabase

  lazy val engineDatabaseInterface: EngineSqlDatabase =
    EngineSlickDatabase.fromParentConfig().initialized(EngineServicesStore.EngineLiquibaseSettings)
}
