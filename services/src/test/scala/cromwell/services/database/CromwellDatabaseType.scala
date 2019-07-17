package cromwell.services.database

import com.typesafe.config.Config
import cromwell.database.migration.liquibase.LiquibaseSettings
import cromwell.database.slick.{EngineSlickDatabase, MetadataSlickDatabase, SlickDatabase}
import cromwell.services.{EngineServicesStore, MetadataServicesStore}

/**
  * The type of database, either Engine or Metadata.
  */
sealed trait CromwellDatabaseType[T <: SlickDatabase] {
  val name: String
  val liquibaseSettings: LiquibaseSettings

  def newDatabase(config: Config): T

  override def toString: String = name
}

object CromwellDatabaseType {
  val All: Seq[CromwellDatabaseType[_ <: SlickDatabase]] = List(
    EngineDatabaseType,
    MetadataDatabaseType,
  )
}

object EngineDatabaseType extends CromwellDatabaseType[EngineSlickDatabase] {
  override val name: String = "Engine"
  override val liquibaseSettings: LiquibaseSettings = EngineServicesStore.EngineLiquibaseSettings

  override def newDatabase(config: Config): EngineSlickDatabase = new EngineSlickDatabase(config)
}

object MetadataDatabaseType extends CromwellDatabaseType[MetadataSlickDatabase] {
  override val name: String = "Metadata"
  override val liquibaseSettings: LiquibaseSettings = MetadataServicesStore.MetadataLiquibaseSettings

  override def newDatabase(config: Config): MetadataSlickDatabase = new MetadataSlickDatabase(config)
}
