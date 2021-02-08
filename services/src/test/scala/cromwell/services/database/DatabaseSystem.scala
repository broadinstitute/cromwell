package cromwell.services.database

/**
  * Cromwell unit tested DBMS. Each DBMS must match a database spun up in test.inc.sh.
  */
// Someday https://github.com/lloydmeta/enumeratum, someday...
sealed trait DatabaseSystem {
  val name: String
  val platform: DatabasePlatform

  override def toString: String = name
}

object DatabaseSystem {
  val All: Seq[DatabaseSystem] = List(
    HsqldbDatabaseSystem,
    MariadbEarliestDatabaseSystem,
    MariadbLatestDatabaseSystem,
    MysqlEarliestDatabaseSystem,
    MysqlLatestDatabaseSystem,
    PostgresqlEarliestDatabaseSystem,
    PostgresqlLatestDatabaseSystem,
  )
}

case object HsqldbDatabaseSystem extends DatabaseSystem {
  override val name: String = "HSQLDB"
  override val platform: HsqldbDatabasePlatform.type = HsqldbDatabasePlatform
}

sealed trait NetworkDatabaseSystem extends DatabaseSystem {
  val dockerImageVersion: String
}

case object MariadbEarliestDatabaseSystem extends NetworkDatabaseSystem {
  override val name: String = "MariaDB"
  override val platform: MariadbDatabasePlatform.type = MariadbDatabasePlatform
  override val dockerImageVersion: String = "5.5"
}

case object MariadbLatestDatabaseSystem extends NetworkDatabaseSystem {
  override val name: String = "MariaDB (latest)"
  override val platform: MariadbDatabasePlatform.type = MariadbDatabasePlatform
  override val dockerImageVersion: String = "latest"
}

case object MysqlEarliestDatabaseSystem extends NetworkDatabaseSystem {
  override val name: String = "MySQL"
  override val platform: MysqlDatabasePlatform.type = MysqlDatabasePlatform
  override val dockerImageVersion: String = "5.6"
}

case object MysqlLatestDatabaseSystem extends NetworkDatabaseSystem {
  override val name: String = "MySQL (latest)"
  override val platform: MysqlDatabasePlatform.type = MysqlDatabasePlatform
  override val dockerImageVersion: String = "latest"
}

case object PostgresqlEarliestDatabaseSystem extends NetworkDatabaseSystem {
  override val name: String = "PostgreSQL"
  override val platform: PostgresqlDatabasePlatform.type = PostgresqlDatabasePlatform
  override val dockerImageVersion: String = "9.5"
}

case object PostgresqlLatestDatabaseSystem extends NetworkDatabaseSystem {
  override val name: String = "PostgreSQL (latest)"
  override val platform: PostgresqlDatabasePlatform.type = PostgresqlDatabasePlatform
  override val dockerImageVersion: String = "latest"
}
