package cromwell.database.slick.tables

import slick.driver.JdbcProfile

trait DriverComponent {
  val driver: JdbcProfile
}
