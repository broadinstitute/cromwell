package cromwell.engine.db.slick

import java.sql.DriverManager

import slick.jdbc.SimpleJdbcAction

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.language.postfixOps

object RunMysql {

  def main(args: Array[String]) = jdbcMain(args)

  def jdbcMain(args: Array[String]): Unit = {
    Class.forName("com.mysql.jdbc.Driver").newInstance()
    val config = DatabaseConfig.databaseConfig

    val Seq(url, user, password) = Seq("url", "user", "password") map {
      key => if (config.hasPath(key)) config.getString(key) else ""
    }

    val connectionString = s"$url?useSSL=true&user=$user&password=$password"
    println(s"Testing connection '$connectionString'")
    val connection = DriverManager.getConnection(connectionString)
    println("Connection valid? " + connection.isValid(5))
    val resultSet = connection.createStatement().executeQuery("select count(*) from WORKFLOW_EXECUTION")
    val first = resultSet.next()
    if (!first) throw new RuntimeException("No first row!")
    val workflowExecutions = resultSet.getInt(1)
    println(s"Found $workflowExecutions workflow executions.")
  }

  def slickMain(args: Array[String]): Unit = {
    val dataAccess = new SlickDataAccess()

    val connectionValidFuture = dataAccess.database.run(SimpleJdbcAction(_.connection.isValid(1)))

    println("Connection valid? " + Await.result(connectionValidFuture, 5 seconds))
  }
}
