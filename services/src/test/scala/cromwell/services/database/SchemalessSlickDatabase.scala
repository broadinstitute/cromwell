package cromwell.services.database

import com.typesafe.config.Config
import cromwell.database.slick.SlickDatabase
import cromwell.database.slick.tables.DataAccessComponent

/**
  * Connects to the passed in config, but does not use nor initialize the schema.
  */
class SchemalessSlickDatabase(config: Config) extends SlickDatabase(config) {
  override lazy val dataAccess: DataAccessComponent = new DataAccessComponent {
    override lazy val driver = slickConfig.profile
    override lazy val schema = driver.DDL("", "")
  }
}
