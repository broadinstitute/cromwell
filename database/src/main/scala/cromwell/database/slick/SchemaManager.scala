package cromwell.database.slick

import java.sql.Connection

import com.typesafe.config.Config
import cromwell.database.liquibase.LiquibaseUtils
import lenthall.config.ScalaConfig._
import liquibase.diff.DiffResult
import slick.dbio.DBIO
import slick.driver.JdbcProfile

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, ExecutionContext, Future}

object SchemaManager {
  /**
    * Returns a new schema manager based on the config key "schema.manager". If the config value is "liquibase", this
    * method returns a new LiquibaseSchemaManager. If the config value is "slick", this method returns a new
    * SlickSchemaManager. If the config key "schema.manager" is not present, returns a new LiquibaseSchemaManager.
    *
    * @param databaseConfig The database configuration information, with a the "schema.manager" key.
    * @return Schema manager.
    */
  def fromConfig(databaseConfig: Config): SchemaManager = {
    databaseConfig.getStringOr("schema.manager", "liquibase").toLowerCase match {
      case "liquibase" => new LiquibaseSchemaManager
      case "slick" => new SlickSchemaManager
    }
  }

  /**
    * Lends a connection to a block of code.
    *
    * @param profile The slick jdbc profile for accessing the database.
    * @param database The database to use for the connection.
    * @param block The block of code to run over the connection.
    * @tparam Profile The slick jdbc profile for accessing the database.
    * @tparam T The return type of the block.
    * @return The return value of the block.
    */
  def withConnection[Profile <: JdbcProfile, T](profile: Profile, database: Profile#Backend#Database)
                                               (block: Connection => T): T = {
    /*
     TODO: Should this withConnection() method have a (implicit?) timeout parameter, that it passes on to Await.result?
     If we run completely asynchronously, nest calls to withConnection, and then call flatMap, the outer connection may
     already be closed before an inner block finishes running.
     */
    Await.result(database.run(profile.api.SimpleDBIO(context => block(context.connection))), Duration.Inf)
  }

  /**
    * Lends two connections to a block of code.
    *
    * @param profile1  The slick jdbc profile for accessing the first database.
    * @param database1 The database to use for the first connection.
    * @param profile2  The slick jdbc profile for accessing the second database.
    * @param database2 The database to use for the second connection.
    * @param block     The block of code to run over the first and second connections.
    * @tparam Profile1 The slick jdbc profile for accessing the first database.
    * @tparam Profile2 The slick jdbc profile for accessing the second database.
    * @tparam T        The return type of the block.
    * @return The return value of the block.
    */
  def withConnections[Profile1 <: JdbcProfile, Profile2 <: JdbcProfile, T]
  (profile1: Profile1, database1: Profile1#Backend#Database, profile2: Profile2, database2: Profile2#Backend#Database)
  (block: (Connection, Connection) => T): T = {
    withConnection(profile1, database1) { connection1 =>
      withConnection(profile2, database2) { connection2 =>
        block(connection1, connection2)
      }
    }
  }
}

sealed trait SchemaManager {
  /**
    * Updates the schema. May require that the schema be empty before running, or may support incremental updates.
    *
    * @param profile The slick jdbc profile for accessing the database.
    * @param schema The schema to create.
    * @tparam Profile The slick jdbc profile for accessing the database.
    * @return An asynchronous runner for the update, possibly eventually returning an error.
    */
  def updateSchema[Profile <: JdbcProfile](profile: Profile, schema: Profile#SchemaDescription): DBIO[Unit]

  /**
    * Updates the schema using the manager specified in the database configuration.
    *
    * May require that the schema be empty before running, or may support incremental updates.
    *
    * @param profile The slick jdbc profile for accessing the database.
    * @param schema The schema to create.
    * @param database The database to update.
    * @tparam Profile The slick jdbc profile for accessing the database.
    * @return An asynchronous runner for the update, possibly eventually returning an error.
    */
  def updateSchema[Profile <: JdbcProfile](profile: Profile,
                                           schema: Profile#SchemaDescription,
                                           database: Profile#Backend#Database): Future[Unit] = {
    database.run(updateSchema(profile, schema))
  }
}

object LiquibaseSchemaManager {

  def compare[ReferenceProfile <: JdbcProfile, ComparisonProfile <: JdbcProfile]
  (referenceProfile: ReferenceProfile,
   referenceDatabase: ReferenceProfile#Backend#Database,
   comparisonProfile: ComparisonProfile,
   comparisonDatabase: ComparisonProfile#Backend#Database)
  (implicit executor: ExecutionContext): DiffResult = {

    import SchemaManager._

    withConnections(referenceProfile, referenceDatabase, comparisonProfile, comparisonDatabase)(LiquibaseUtils.compare)
  }
}

/**
  * Creates or updates a schema using liquibase.
  */
class LiquibaseSchemaManager extends SchemaManager {
  override def updateSchema[Profile <: JdbcProfile](profile: Profile, schema: Profile#SchemaDescription) = {
    import profile.api._
    SimpleDBIO(context => LiquibaseUtils.updateSchema(context.connection))
  }
}

/**
  * Creates a schema using Slick's built in utility.
  *
  * NOTE: The schema must not be already present or the update will fail.
  */
class SlickSchemaManager extends SchemaManager {
  override def updateSchema[Profile <: JdbcProfile](profile: Profile, schema: Profile#SchemaDescription) = {
    profile.createSchemaActionExtensionMethods(schema.asInstanceOf[profile.SchemaDescription]).create
  }
}
