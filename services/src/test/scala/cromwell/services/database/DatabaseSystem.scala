package cromwell.services.database

/**
  * Cromwell supported DBMS.
  */
// Someday https://github.com/lloydmeta/enumeratum, someday...
sealed trait DatabaseSystem {
  val productName: String
  val shortName: String
  val configPath: String

  override def toString: String = productName
}

object DatabaseSystem {
  def apply(productName: String): DatabaseSystem = {
    productName match {
      case MysqlDatabaseSystem.productName => MysqlDatabaseSystem
      case HsqldbDatabaseSystem.productName => HsqldbDatabaseSystem
      case PostgresqlDatabaseSystem.productName => PostgresqlDatabaseSystem
      case MariadbDatabaseSystem.productName => MariadbDatabaseSystem
      case _ => throw new UnsupportedOperationException(s"Unknown database system: $productName")
    }
  }

  val All: Seq[DatabaseSystem] = List(
    HsqldbDatabaseSystem,
    MariadbDatabaseSystem,
    MysqlDatabaseSystem,
    PostgresqlDatabaseSystem,
  )
}

case object HsqldbDatabaseSystem extends DatabaseSystem {
  override val productName: String = "HSQL Database Engine"
  override val shortName: String = "HSQLDB"
  override val configPath: String = "database"
}

case object MariadbDatabaseSystem extends DatabaseSystem {
  override val productName: String = "MariaDB"
  override val shortName = productName
  override val configPath: String = "database-test-mariadb"
}

case object MysqlDatabaseSystem extends DatabaseSystem {
  override val productName: String = "MySQL"
  override val shortName = productName
  override val configPath: String = "database-test-mysql"
}

case object PostgresqlDatabaseSystem extends DatabaseSystem {
  override val productName: String = "PostgreSQL"
  override val shortName = productName
  override val configPath: String = "database-test-postgresql"
}
