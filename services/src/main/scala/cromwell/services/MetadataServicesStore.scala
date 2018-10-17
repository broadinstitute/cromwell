package cromwell.services

import cromwell.database.migration.liquibase.LiquibaseSettings
import cromwell.database.slick.MetadataSlickDatabase
import cromwell.database.sql.{MetadataSqlDatabase, StreamMetadataSqlDatabase}
import liquibase.database.Database

trait MetadataServicesStore {
  def metadataDatabaseInterface: MetadataSqlDatabase with StreamMetadataSqlDatabase = MetadataServicesStore.metadataDatabaseInterface
}

object MetadataServicesStore {
  val MetadataLiquibaseSettings = LiquibaseSettings(
    "sql_metadata_changelog.xml",
    "SQLMETADATA" + Database.databaseChangeLogLockTableName,
    "SQLMETADATA" + Database.databaseChangeLogTableName
  )

  import ServicesStore.EnhancedSqlDatabase

  // Mix in AutoCloseable for shutdown purposes.
  // This whole MetadataSqlDatabase interface will likely change to a more abstract MetadataDAO in the future.
  val metadataDatabaseInterface: MetadataSqlDatabase with StreamMetadataSqlDatabase =
    MetadataSlickDatabase.fromParentConfig().initialized(MetadataServicesStore.MetadataLiquibaseSettings)
}
