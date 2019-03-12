package cromwell.database.slick

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.database.slick.tables.EngineDataAccessComponent
import cromwell.database.sql.EngineSqlDatabase

import scala.concurrent.duration.Duration

object EngineSlickDatabase {
  def fromParentConfig(parentConfig: Config = ConfigFactory.load): EngineSlickDatabase = {
    val (databaseConfig, queryTimeout) = SlickDatabase.getDatabaseConfig("engine", parentConfig)
    new EngineSlickDatabase(databaseConfig, queryTimeout)
  }
}

class EngineSlickDatabase(originalDatabaseConfig: Config, queryTimeout: Duration)
  extends SlickDatabase(originalDatabaseConfig, queryTimeout)
    with EngineSqlDatabase
    with WorkflowStoreSlickDatabase
    with JobKeyValueSlickDatabase
    with JobStoreSlickDatabase
    with CallCachingSlickDatabase
    with SubWorkflowStoreSlickDatabase
    with DockerHashStoreSlickDatabase {
  override lazy val dataAccess = new EngineDataAccessComponent(slickConfig.profile)
}
