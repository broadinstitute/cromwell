package cromwell.engine.db.slick

import com.typesafe.config.Config
import cromwell.engine.db.LiquibaseUtils
import slick.dbio.DBIO
import slick.driver.JdbcProfile
import lenthall.config.ScalaConfig._

object SchemaManager {
  /**
    * Updates the schema using the manager specified in the database configuration.
    *
    * The configuration key is "schemaManager", and may have a value of "liquibase" or "slick". If left blank uses
    * "liquibase".
    *
    * May require that the schema be empty before running, or may support incremental updates.
    *
    * @param databaseConfig The database configuration information, including the "schemaManager" key.
    * @param profile The slick jdbc profile for accessing the database.
    * @param schema The schema to create.
    * @tparam Profile The slick jdbc profile for accessing the database.
    * @return An asynchronous runner for the update, possibly eventually returning an error.
    */
  def updateSchema[Profile <: JdbcProfile](databaseConfig: Config,
                                           profile: Profile,
                                           schema: Profile#SchemaDescription): DBIO[Unit] = {
    val schemaManager = databaseConfig.getStringOr("schema.manager", "liquibase").toLowerCase match {
      case "liquibase" => new LiquibaseSchemaManager
      case "slick" => new SlickSchemaManager
    }
    schemaManager.updateSchema(profile, schema)
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
