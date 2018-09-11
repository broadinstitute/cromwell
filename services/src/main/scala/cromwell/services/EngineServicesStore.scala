package cromwell.services

import cromwell.database.migration.liquibase.LiquibaseSettings
import cromwell.database.slick.EngineSlickDatabase
import cromwell.database.sql.EngineSqlDatabase

object EngineServicesStore {
  val EngineLiquibaseSettings = LiquibaseSettings("changelog.xml")

  import ServicesStore.EnhancedSqlDatabase

  val engineDatabaseInterface: EngineSqlDatabase =
    EngineSlickDatabase.fromParentConfig().initialized(EngineServicesStore.EngineLiquibaseSettings)
}
