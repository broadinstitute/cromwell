package cromwell.engine.db.slick

import java.sql.Connection

import cromwell.util.DatabaseConfig
import liquibase.Liquibase
import liquibase.database.DatabaseConnection
import liquibase.resource.{FileSystemResourceAccessor, ResourceAccessor}

// Modified from https://tillias.wordpress.com/2012/11/10/unit-testing-and-integration-testing-using-junit-liquibase-hsqldb-hibernate-maven-and-spring-framework/
trait SlickTestDatabase {
  SlickTestDatabase.checkInitialized()
}

object SlickTestDatabase {
  def checkInitialized() {
    // do nothing, static constructor run of start() does actual work
  }

  start()

  private def start(): Unit = {
    if (DatabaseConfig.liquibaseSetup)
      setUp("test")
  }

  private def setUp(contexts: String) {
    val connectionClass = Class.forName(DatabaseConfig.liquibaseConnection)
    val connectionConstructor = connectionClass.getConstructor(classOf[Connection])
    val jdbcConnection = DataAccessController.database.source.createConnection()
    val liquibaseConnection: DatabaseConnection =
      connectionConstructor.newInstance(jdbcConnection).asInstanceOf[DatabaseConnection]
    try {
      val resourceAccessor: ResourceAccessor = new FileSystemResourceAccessor()
      val liquibase = new Liquibase(DatabaseConfig.liquibaseChangeLog, resourceAccessor, liquibaseConnection)
      if (DatabaseConfig.liquibaseDropAll)
        liquibase.dropAll()
      liquibase.update(contexts)
    } finally {
      liquibaseConnection.close()
    }
  }
}
