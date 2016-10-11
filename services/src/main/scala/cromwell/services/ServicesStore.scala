package cromwell.services

import com.typesafe.config.ConfigFactory
import cromwell.database.migration.liquibase.LiquibaseUtils
import cromwell.database.slick.SlickDatabase
import cromwell.database.sql.SqlDatabase
import net.ceedubs.ficus.Ficus._
import org.slf4j.LoggerFactory

trait ServicesStore {
  def databaseInterface: SqlDatabase
}

object ServicesStore {

  implicit class EnhancedSqlDatabase[A <: SqlDatabase](val sqlDatabase: A) extends AnyVal {
    def initialized: A = {
      if (sqlDatabase.databaseConfig.as[Option[Boolean]]("liquibase.updateSchema").getOrElse(true)) {
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

  private lazy val log = LoggerFactory.getLogger("SingletonServicesStore")
  private val databaseConfig = {
    val config = ConfigFactory.load.getConfig("database")
    if (config.hasPath("config")) {
      log.warn(
        """
          |Use of configuration path 'database.config' is deprecated.
          |
          |Move the configuration directly under the 'database' element, and remove the key 'database.config'.
          |""".stripMargin)
      config.getConfig(config.getString("config"))
    } else {
      config
    }
  }

  import ServicesStore.EnhancedSqlDatabase

  val databaseInterface: SqlDatabase = new SlickDatabase(databaseConfig).initialized
}
