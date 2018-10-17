package centaur

import com.typesafe.config.Config
import cromwell.database.slick.{EngineSlickDatabase, MetadataSlickDatabase}
import cromwell.database.sql.{EngineSqlDatabase, MetadataSqlDatabase, StreamMetadataSqlDatabase}

trait CromwellDatabase {
  def engineDatabase: EngineSqlDatabase
  def metadataDatabase: MetadataSqlDatabase
}

object CromwellDatabase {
  /**
    * Wraps connections to a cromwell database. The database connections are not initialized until first use.
    */
  def fromConfig(rootConfig: Config): CromwellDatabase = {
    lazy val cromwellConfig = rootConfig.getConfig("cromwell")
    new CromwellDatabase {
      override lazy val engineDatabase: EngineSqlDatabase = EngineSlickDatabase.fromParentConfig(cromwellConfig)
      override lazy val metadataDatabase: MetadataSqlDatabase with StreamMetadataSqlDatabase = MetadataSlickDatabase.fromParentConfig(cromwellConfig)
    }
  }
}
