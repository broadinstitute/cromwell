package centaur

import java.net.URL
import java.nio.file.{Path, Paths}
import java.util.concurrent.TimeUnit

import com.typesafe.config.{Config, ConfigFactory}
import configs.syntax._

import scala.concurrent.duration.FiniteDuration

object CentaurRunMode {
  def apply(config: Config): CentaurRunMode = {
    val cromwellConfig = config.getConfig("cromwell")
    cromwellConfig.get[String]("mode").value match {
      case "url" =>
        val url = cromwellConfig.get[String]("url").value
        UnmanagedCromwellServer(new URL(url))
      case "jar" =>
        val jarConf = cromwellConfig.get[Config]("jar").value
        val preRestart = JarCromwellConfiguration(jarConf)
        val withRestart = jarConf.getBoolean("withRestart")
        val postRestartConfig =
          JarCromwellConfiguration(cromwellConfig.get[Config]("post-restart-jar").valueOrElse(jarConf))

        ManagedCromwellServer(preRestart, postRestartConfig, withRestart)
      case "docker-compose" =>
        val composeConf = cromwellConfig.get[Config](path = "docker-compose").value
        val preRestart = DockerComposeCromwellConfiguration(composeConf)
        val withRestart = composeConf.getBoolean("withRestart")
        val postRestartConfig =
          DockerComposeCromwellConfiguration(cromwellConfig.getOrElse("post-restart-docker-compose", composeConf).value)
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

object CentaurConfig {
  lazy val conf: Config = ConfigFactory.load().getConfig("centaur")
  
  lazy val runMode: CentaurRunMode = CentaurRunMode(conf)
  lazy val expectCarbonite: Boolean = conf.getOrElse("expectCarbonite", false).value
  
  lazy val cromwellUrl: URL = runMode.cromwellUrl
  lazy val sendReceiveTimeout: FiniteDuration = conf.getDuration("sendReceiveTimeout").toScala
  lazy val maxWorkflowLength: FiniteDuration = conf.getDuration("maxWorkflowLength").toScala
  lazy val metadataConsistencyTimeout: FiniteDuration = conf.getDuration("metadataConsistencyTimeout").toScala

  lazy val metadataDeletionMinimumWait: FiniteDuration = conf.getDuration("metadataDeletionMinimumWait").toScala
  lazy val metadataDeletionMaximumWait: FiniteDuration = conf.getDuration("metadataDeletionMaximumWait").toScala

  lazy val standardTestCasePath: Path = Paths.get(conf.getString("standardTestCasePath"))

  // If provided, any tests will be appended to the tests in standardTestCasePath
  lazy val optionalTestPath: Option[Path] = conf.get[Option[Path]]("optionalTestPath").value

  implicit class EnhancedJavaDuration(val javaDuration: java.time.Duration) extends AnyVal {
    def toScala: FiniteDuration = FiniteDuration(javaDuration.toMillis, TimeUnit.MILLISECONDS)
  }
}
