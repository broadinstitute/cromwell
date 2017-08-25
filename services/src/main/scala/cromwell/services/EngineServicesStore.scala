package cromwell.services

import cromwell.database.migration.liquibase.LiquibaseSettings
import cromwell.database.slick.EngineSlickDatabase
import cromwell.database.sql.EngineSqlDatabase

object EngineServicesStore {
  val EngineLiquibaseSettings = LiquibaseSettings("changelog.xml")

  private val engineDatabaseConfig = ServicesStore.getDatabaseConfig("engine")

  import ServicesStore.EnhancedSqlDatabase

  val engineDatabaseInterface: EngineSqlDatabase =
    new EngineSlickDatabase(engineDatabaseConfig)
      .initialized(EngineServicesStore.EngineLiquibaseSettings)
}
