package cromwell.engine.db.slick

import java.sql.{Connection, DriverManager}

import com.typesafe.config.{ConfigFactory, ConfigRenderOptions}
import slick.jdbc.SimpleJdbcAction

import scala.collection.JavaConverters._
import scala.concurrent.Await
import scala.concurrent.duration._
import scala.language.postfixOps

object RunMysql {
  Class.forName("com.mysql.jdbc.Driver").newInstance()

  val user: String = getConfig("user")
  val password: String = getConfig("password")
  val url: String = getConfig("url")
  val urlWithSsl: String = withSSL(url)
  val urlWithoutSsl: String = stripSSL(url)
  val configWithoutConnection = DatabaseConfig.databaseConfig.
    withoutPath("url").
    withoutPath("driver").
    withoutPath("dataSourceClass").
    withoutPath("properties.url").
    withoutPath("properties.useSSL").
    withoutPath("properties.requireSSL")

  var successes: Seq[String] = Seq.empty
  var failures: Seq[String] = Seq.empty

  // Hiding this main() from sbt, but leaving class for future debugging.
  //def main(args: Array[String]) = runMain()

  def runMain() = {
    println("user: " + user)
    println("password: " + ("*" * password.length))
    println("url: " + url)
    println("urlWithSsl: " + urlWithSsl)
    println("urlWithoutSsl: " + urlWithoutSsl)

    tryRun("jdbcMain", jdbcMain())
    tryRun("slickMain", slickMain())
    tryRun("jdbcRequireSsl", jdbcRequireSSL())

    tryRun("dataSourceRaw", dataSourceTest(urlWithoutSsl, Option(false)))
    tryRun("dataSourceRawProp", dataSourceTest(urlWithoutSsl, Option(true)))
    tryRun("dataSourceRawNoProp", dataSourceTest(urlWithoutSsl, None))

    tryRun("dataSourceSslUrl", dataSourceTest(urlWithSsl, Option(false)))
    tryRun("dataSourceSslUrlProp", dataSourceTest(urlWithSsl, Option(true)))
    tryRun("dataSourceSslNoProp", dataSourceTest(urlWithSsl, None))

    tryRun("slickRaw", slickTest(urlWithoutSsl, sslProp = Option(false)))
    tryRun("slickRawProp", slickTest(urlWithoutSsl, sslProp = Option(true)))
    tryRun("slickRawNoProp", slickTest(urlWithoutSsl, sslProp = None))
    tryRun("slickRawDriver", slickTest(urlWithoutSsl, dataSource = false))

    tryRun("slickSslUrl", slickTest(urlWithSsl, sslProp = Option(false)))
    tryRun("slickSslUrlProp", slickTest(urlWithSsl, sslProp = Option(true)))
    tryRun("slickSslUrlNoProp", slickTest(urlWithSsl, sslProp = None))
    tryRun("slickSslDriver", slickTest(urlWithSsl, dataSource = false))

    println("Successes: " + successes.size)
    successes.foreach(s => println("  " + s))
    println("Failures: " + failures.size)
    failures.foreach(f => println("  " + f))
  }

  // Tests

  def jdbcMain() = {
    val connectionString = s"$url?useSSL=true"
    println(s"Testing connection '$connectionString'")
    val connection = DriverManager.getConnection(connectionString, user, password)
    testConnection(connection)
  }

  def jdbcRequireSSL() = {
    val connectionString = s"$url?useSSL=true&requireSSL=true"
    println(s"Testing connection '$connectionString'")
    val connection = DriverManager.getConnection(connectionString, user, password)
    testConnection(connection)
  }

  def slickMain(): Unit = {
    val dataAccess = new SlickDataAccess()
    val connectionFuture = dataAccess.database.run(SimpleJdbcAction(context => testConnection(context.connection)))
    Await.result(connectionFuture, 10 seconds)
  }

  def dataSourceMain(): Unit = {
    val ds = new com.mysql.jdbc.jdbc2.optional.MysqlDataSource
    ds.setURL(url)
    ds.setUseSSL(true)
    ds.setRequireSSL(true)
    testConnection(ds.getConnection(user, password))
  }

  // Datasource permutations

  def dataSourceTest(testUrl: String, useSsl: Option[Boolean]): Unit = {
    val ds = new com.mysql.jdbc.jdbc2.optional.MysqlDataSource
    ds.setURL(testUrl)
    useSsl foreach ds.setUseSSL
    useSsl foreach ds.setRequireSSL
    testConnection(ds.getConnection(user, password))
  }

  // Slick tests

  def slickTest(urlVal: String,
            dataSource: Boolean = true,
            sslProp: Option[Boolean] = None): Unit = {

    var map: Map[String, AnyRef] = Map.empty

    if (dataSource) {
      map += "dataSourceClass" -> "com.mysql.jdbc.jdbc2.optional.MysqlDataSource"
      map += "properties.url" -> urlVal
      sslProp foreach { useSsl =>
        map += "properties.useSSL" -> Boolean.box(useSsl)
        map += "properties.requireSSL" -> Boolean.box(useSsl)
      }
    } else {
      map += "driver" -> "com.mysql.jdbc.Driver"
      map += "url" -> urlVal
    }

    val config = ConfigFactory.parseMap(map.asJava, "debug slick options").withFallback(configWithoutConnection)
    println("config:")
    println(config.root().render(ConfigRenderOptions.defaults()))
    val dataAccess = new SlickDataAccess(config)
    val connectionFuture = dataAccess.database.run(SimpleJdbcAction(context => testConnection(context.connection)))
    Await.result(connectionFuture, 10 seconds)
  }

  // Utilities

  private def getConfig(key: String) =
    getConfigOpt("properties." + key) orElse getConfigOpt(key) getOrElse ""

  private def getConfigOpt(key: String): Option[String] = {
    if (DatabaseConfig.databaseConfig.hasPath(key)) {
      Option(DatabaseConfig.databaseConfig.getString(key))
    } else {
      None
    }
  }

  private def tryRun(label: String, test: => Unit): Unit = {
    try {
      println("Running: " + label)
      test
      println("Success: " + label)
      successes :+= label
    } catch {
      case e: Exception =>
        failures :+= label
        println("Failed: " + label)
        e.printStackTrace()
    }
    println()
  }

  private def testConnection(connection: Connection) = {
    try {
      println("Connection url: " + connection.getMetaData.getURL)
      println("Connection valid? " + connection.isValid(5))
      val resultSet = connection.createStatement().executeQuery("select count(*) from WORKFLOW_EXECUTION")
      val first = resultSet.next()
      if (!first) throw new RuntimeException("No first row!")
      val workflowExecutions = resultSet.getInt(1)
      println(s"Found $workflowExecutions workflow executions.")
    } finally{
      try {
        connection.close()
      } catch {
        case e: Exception =>
        /* ignore */
      }
    }
  }

  private def stripSSL(url: String): String = {
    var newUrl = url
    newUrl = newUrl.replaceAll( """\?useSSL=true\&""", "?")
    newUrl = newUrl.replaceAll( """\?useSSL=true$""", "")
    newUrl = newUrl.replaceAll( """\&useSSL=true""", "")
    newUrl = newUrl.replaceAll( """\?requireSSL=true\&""", "?")
    newUrl = newUrl.replaceAll( """\?requireSSL=true$""", "")
    newUrl = newUrl.replaceAll( """\&requireSSL=true""", "")
    newUrl
  }

  private def withSSL(url: String): String = {
    val newUrl = stripSSL(url)
    if (newUrl.contains("?"))
      newUrl + "&useSSL=true&requireSSL=true"
    else
      newUrl + "?useSSL=true&requireSSL=true"
  }
}
