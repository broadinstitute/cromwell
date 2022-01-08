package cromwell.services.database

/**
  * Schema management by either Liquibase (the default) or Slick.
  */
sealed trait SchemaManager {
  val name: String
  def other: SchemaManager

  override def toString: String = name
}

object SchemaManager {
  val All: Seq[SchemaManager] = List(
    LiquibaseSchemaManager,
    SlickSchemaManager,
  )
}

case object LiquibaseSchemaManager extends SchemaManager {
  override val name: String = "Liquibase"
  lazy val other: SchemaManager = SlickSchemaManager
}

case object SlickSchemaManager extends SchemaManager {
  override val name: String = "Slick"
  lazy val other: SchemaManager = LiquibaseSchemaManager
}
