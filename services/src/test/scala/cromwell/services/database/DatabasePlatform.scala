package cromwell.services.database

/**
  * Cromwell supported DBMS platforms.
  */
sealed trait DatabasePlatform {
  def name: String

  override def toString: String = name
}

case object HsqldbDatabasePlatform extends DatabasePlatform {
  override val name: String = "HSQLDB"
}

case object MariadbDatabasePlatform extends DatabasePlatform {
  override val name: String = "MariaDB"
}

case object MysqlDatabasePlatform extends DatabasePlatform {
  override val name: String = "MySQL"
}

case object PostgresqlDatabasePlatform extends DatabasePlatform {
  override val name: String = "PostgreSQL"
}
