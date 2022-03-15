package cromwell.engine.workflow

import com.typesafe.config.Config
import cromwell.CromwellTestKitSpec
import cromwell.database.slick.{EngineSlickDatabase, MetadataSlickDatabase}
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore
import cromwell.services.{EngineServicesStore, MetadataServicesStore}
import cromwell.services.ServicesStore.EnhancedSqlDatabase

import scala.util.Try

trait SqlWorkflowStoreBuilder {

  val rootConfig = CromwellTestKitSpec.DefaultConfig
  val databaseConfig = rootConfig.getConfig("database")

  def runWithDatabase[T](databaseConfig: Config)(block: SqlWorkflowStore => T): T = {
    val database = new EngineSlickDatabase(databaseConfig).initialized(EngineServicesStore.EngineLiquibaseSettings)
    val metadataDatabase = new MetadataSlickDatabase(databaseConfig).initialized(MetadataServicesStore.MetadataLiquibaseSettings)
    try {
      block(SqlWorkflowStore(database, metadataDatabase))
    } finally {
      Try(database.close())
      ()
    }
  }
}
