package cromwell.database.slick

import slick.driver.JdbcProfile

trait DriverComponent {
  val driver: JdbcProfile
}
