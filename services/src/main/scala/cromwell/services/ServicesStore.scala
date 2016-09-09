package cromwell.services

import cromwell.database.core.SqlConfiguration
import cromwell.database.migration.liquibase.LiquibaseUtils
import cromwell.database.slick.SlickDatabase
import cromwell.database.sql.SqlDatabase
import lenthall.config.ScalaConfig._

trait ServicesStore {
  def databaseInterface: SqlDatabase
}

object ServicesStore {

  implicit class EnhancedSqlDatabase[A <: SqlDatabase](val sqlDatabase: A) extends AnyVal {
    def initialized: A = {
      if (sqlDatabase.databaseConfig.getBooleanOr("liquibase.updateSchema", default = true)) {
        sqlDatabase withConnection LiquibaseUtils.updateSchema
      }
      sqlDatabase
    }
  }

}

trait SingletonServicesStore extends ServicesStore {
  final override val databaseInterface: SqlDatabase = SingletonServicesStore.databaseInterface
}

object SingletonServicesStore {

  import ServicesStore.EnhancedSqlDatabase

  val databaseInterface: SqlDatabase = new SlickDatabase(SqlConfiguration.defaultDatabaseConfig).initialized
}
