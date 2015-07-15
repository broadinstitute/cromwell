package cromwell.engine.db.slick

import slick.jdbc.SimpleJdbcAction

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.language.postfixOps

object RunMysql {

  def main(args: Array[String]): Unit = {
    val dataAccess = new SlickDataAccess()

    val connectionValidFuture = dataAccess.database.run(SimpleJdbcAction(_.connection.isValid(1)))

    println("Connection valid? " + Await.result(connectionValidFuture, 5 seconds))
  }
}
