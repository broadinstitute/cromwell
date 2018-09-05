package cromwell.services

import cromwell.database.migration.liquibase.{LiquibaseSettings, LiquibaseUtils}
import cromwell.database.sql.SqlDatabase
import net.ceedubs.ficus.Ficus._

object ServicesStore {
  implicit class EnhancedSqlDatabase[A <: SqlDatabase](val sqlDatabase: A) extends AnyVal {
    def initialized(settings: LiquibaseSettings): A = {
      if (sqlDatabase.databaseConfig.getOrElse("liquibase.updateSchema", true)) {
        sqlDatabase withConnection LiquibaseUtils.updateSchema(settings)
      }
      sqlDatabase
    }
  }
}
