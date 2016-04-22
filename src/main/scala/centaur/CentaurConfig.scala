package centaur

import java.net.URL
import java.nio.file.{Path, Paths}
import java.util.concurrent.TimeUnit

import com.typesafe.config.ConfigFactory
import lenthall.config.ScalaConfig._

import scala.concurrent.duration.FiniteDuration

object CentaurConfig {
  lazy val conf = ConfigFactory.load()
  lazy val cromwellUrl = new URL(conf.getString("centaur.cromwellUrl"))
  lazy val sendReceiveTimeout = conf.getDuration("centaur.sendReceiveTimeout").toScala
  lazy val maxWorkflowLength = conf.getDuration("centaur.maxWorkflowLength").toScala
  lazy val successfulTestCasePath = Paths.get(conf.getString("centaur.successfulTestCasePath"))
  lazy val failingTestCasePath = Paths.get(conf.getString("centaur.failingTestCasePath"))
  lazy val submissionFailureTestCasePath = Paths.get(conf.getString("centaur.submissionFailureTestCasePath"))

  // If provided, any tests will be appended to the tests in successfulTestCasePath
  lazy val optionalTestPath: Option[Path] = conf.getStringOption("centaur.optionalTestPath") map { Paths.get(_) }

  implicit class EnhancedJavaDuration(val javaDuration: java.time.Duration) extends AnyVal {
    def toScala: FiniteDuration = FiniteDuration(javaDuration.toMillis, TimeUnit.MILLISECONDS)
  }
}
