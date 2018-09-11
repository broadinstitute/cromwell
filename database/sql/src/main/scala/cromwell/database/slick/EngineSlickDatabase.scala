package cromwell.database.slick

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.database.slick.tables.EngineDataAccessComponent
import cromwell.database.sql.EngineSqlDatabase

object EngineSlickDatabase {
  def fromParentConfig(parentConfig: Config = ConfigFactory.load): EngineSlickDatabase = {
    val databaseConfig = SlickDatabase.getDatabaseConfig("engine", parentConfig)
    new EngineSlickDatabase(databaseConfig)
  }
}

class EngineSlickDatabase(originalDatabaseConfig: Config)
  extends SlickDatabase(originalDatabaseConfig)
    with EngineSqlDatabase
    with WorkflowStoreSlickDatabase
    with JobKeyValueSlickDatabase
    with JobStoreSlickDatabase
    with CallCachingSlickDatabase
    with SubWorkflowStoreSlickDatabase
    with DockerHashStoreSlickDatabase {
  override lazy val dataAccess = new EngineDataAccessComponent(slickConfig.profile)
}
