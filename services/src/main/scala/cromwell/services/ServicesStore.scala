package cromwell.services

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.database.migration.liquibase.{LiquibaseSettings, LiquibaseUtils}
import cromwell.database.sql.SqlDatabase
import net.ceedubs.ficus.Ficus._

object ServicesStore {

  def getDatabaseConfig(name: String): Config = {
    val rootDatabaseConfig = ConfigFactory.load.getConfig("database")
    val databaseConfig = rootDatabaseConfig.getOrElse(name, rootDatabaseConfig)
    if (databaseConfig.hasPath("config")) {
      val msg =
        """|
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
        """|
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
        new RuntimeException(msg)
    }
    databaseConfig
  }

  implicit class EnhancedSqlDatabase[A <: SqlDatabase](val sqlDatabase: A) extends AnyVal {
    def initialized(settings: LiquibaseSettings): A = {
      if (sqlDatabase.databaseConfig.as[Option[Boolean]]("liquibase.updateSchema").getOrElse(true)) {
        sqlDatabase withConnection LiquibaseUtils.updateSchema(settings)
      }
      sqlDatabase
    }
  }
}
