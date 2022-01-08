package cromwell.services.database

/**
  * Metadata from the JDBC connection and driver.
  */
case class ConnectionMetadata
(
  databaseProductName: String,
  databaseProductVersion: String,
  databaseMajorVersion: Int,
  databaseMinorVersion: Int,
  driverName: String,
  driverVersion: String,
  driverMajorVersion: Int,
  driverMinorVersion: Int,
  jdbcMajorVersion: Int,
  jdbcMinorVersion: Int,
)
