package cromwell.services

import com.typesafe.config.ConfigFactory
import cromwell.database.migration.liquibase.LiquibaseUtils
import cromwell.database.slick.SlickDatabase
import cromwell.database.sql.SqlDatabase
import net.ceedubs.ficus.Ficus._

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

  private val databaseConfig = ConfigFactory.load.getConfig("database")

  if (databaseConfig.hasPath("config")) {
      val msg = """
        |*******************************
        |***** DEPRECATION MESSAGE *****
        |*******************************
        |
        |Use of configuration path 'database.config' has been deprecated.
        |
        |Move the configuration directly under the 'database' element, and remove the key 'database.config'.
        |
        |""".stripMargin
    throw new Exception(msg)
  } else if (databaseConfig.hasPath("driver")) {
    val msg =
      """
      |*******************************
      |***** DEPRECATION MESSAGE *****
      |*******************************
      |
      |Use of configuration path 'database.driver' has been deprecated. Replace with a "profile" element instead, e.g:
      |
      |database {
      |  #driver = "slick.driver.MySQLDriver$" #old
      |  profile = "slick.jdbc.MySQLProfile$"  #new
      |  db {
      |    driver = "com.mysql.jdbc.Driver"
      |    url = "jdbc:mysql://host/cromwell?rewriteBatchedStatements=true"
      |    user = "user"
      |    password = "pass"
      |    connectionTimeout = 5000
      |  }
      |}
      |
      |Cromwell thanks you.
      |""".stripMargin
    throw
      new Exception(msg)
  }

  import ServicesStore.EnhancedSqlDatabase

  val databaseInterface: SqlDatabase = new SlickDatabase(databaseConfig).initialized
}
