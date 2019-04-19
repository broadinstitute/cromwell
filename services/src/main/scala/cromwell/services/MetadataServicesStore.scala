package cromwell.services

import cromwell.database.migration.liquibase.LiquibaseSettings
import cromwell.database.slick.MetadataSlickDatabase
import cromwell.database.sql.MetadataSqlDatabase
import liquibase.database.Database

trait MetadataServicesStore {
  def metadataDatabaseInterface: MetadataSqlDatabase = MetadataServicesStore.metadataDatabaseInterface
}

object MetadataServicesStore {
  lazy val MetadataLiquibaseSettings = LiquibaseSettings(
    "sql_metadata_changelog.xml",
    "SQLMETADATA" + Database.databaseChangeLogLockTableName,
    "SQLMETADATA" + Database.databaseChangeLogTableName
  )

  import ServicesStore.EnhancedSqlDatabase

  // Mix in AutoCloseable for shutdown purposes.
  // This whole MetadataSqlDatabase interface will likely change to a more abstract MetadataDAO in the future.
  lazy val metadataDatabaseInterface: MetadataSqlDatabase =
    MetadataSlickDatabase.fromParentConfig().initialized(MetadataServicesStore.MetadataLiquibaseSettings)
}
