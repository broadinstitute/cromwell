package centaur

import java.net.URL
import java.nio.file.{Path, Paths}
import java.util.concurrent.TimeUnit
import better.files._

import com.typesafe.config.{Config, ConfigFactory}
import configs.syntax._

import scala.concurrent.duration.FiniteDuration

object CromwellConfiguration {
  def apply(conf: Config): CromwellConfiguration = {
    val jarPath = conf.getString("path")
    val confPath = conf.getString("conf")
    val logPath = conf.getString("log")
    
    new CromwellConfiguration(jarPath, confPath, logPath)
  }
}

object CentaurRunMode {
  def apply(config: Config) = {
    val cromwellConfig = config.getConfig("cromwell")
    cromwellConfig.get[String]("mode").value match {
      case "url" =>
        val url = cromwellConfig.get[String]("url").value
        UnmanagedCromwellServer(new URL(url))
      case "jar" =>
        val jarConf = cromwellConfig.get[Config]("jar").value
        val preRestart = CromwellConfiguration(jarConf)
        val withRestart = jarConf.getBoolean("withRestart")
        val postRestartConfig =
          CromwellConfiguration(cromwellConfig.get[Config]("post-restart-jar").valueOrElse(jarConf))

        ManagedCromwellServer(preRestart, postRestartConfig, withRestart)
      case other => throw new Exception(s"Unrecognized cromwell mode: $other")
    }
  }
}

sealed trait CentaurRunMode {
  def cromwellUrl: URL
}

case class UnmanagedCromwellServer(cromwellUrl : URL) extends CentaurRunMode
case class ManagedCromwellServer(preRestart: CromwellConfiguration, postRestart: CromwellConfiguration, withRestart: Boolean) extends CentaurRunMode {
  override val cromwellUrl = new URL(s"http://localhost:${CromwellManager.ManagedCromwellPort}")
}

case class CromwellConfiguration(jar: String, conf: String, logFile: String)

object CentaurConfig {
  lazy val conf = ConfigFactory.load().getConfig("centaur")
  
  lazy val runMode = CentaurRunMode(conf)
  
  lazy val cromwellUrl = runMode.cromwellUrl
  lazy val sendReceiveTimeout = conf.getDuration("sendReceiveTimeout").toScala
  lazy val maxWorkflowLength = conf.getDuration("maxWorkflowLength").toScala
  lazy val metadataConsistencyTimeout = conf.getDuration("metadataConsistencyTimeout").toScala

  lazy val standardTestCasePath = Paths.get(conf.getString("standardTestCasePath"))

  // If provided, any tests will be appended to the tests in standardTestCasePath
  lazy val optionalTestPath: Option[Path] = conf.get[Option[Path]]("optionalTestPath") valueOrElse None
  // If provided, the token will become the default value for the workflow option "refresh_token"
  lazy val optionalTokenPath: Option[Path] = conf.get[Option[Path]]("optionalTokenPath") valueOrElse None
  lazy val optionalToken: Option[String] = optionalTokenPath.map(File(_).contentAsString.trim)

  implicit class EnhancedJavaDuration(val javaDuration: java.time.Duration) extends AnyVal {
    def toScala: FiniteDuration = FiniteDuration(javaDuration.toMillis, TimeUnit.MILLISECONDS)
  }
}
