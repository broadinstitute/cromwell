package centaur

import centaur.CromwellManager.ManagedCromwellPort
import com.typesafe.config.Config

object JarCromwellConfiguration {
  def apply(conf: Config): CromwellConfiguration = {
    val jarPath = conf.getString("path")
    val confPath = conf.getString("conf")
    val logPath = conf.getString("log")

    JarCromwellConfiguration(jarPath, confPath, logPath)
  }
}

case class JarCromwellConfiguration(jar: String, conf: String, logFile: String) extends CromwellConfiguration {
  override def createProcess: CromwellProcess = {
    case class JarCromwellProcess(override val cromwellConfiguration: JarCromwellConfiguration) extends CromwellProcess {
      private val command = Array(
        "java",
        s"-Dconfig.file=$conf",
        s"-Dwebservice.port=$ManagedCromwellPort",
        "-jar",
        jar,
        "server")

      private var process: Option[Process] = None

      override def start(): Unit = {
        process = Option(runProcess(command))
      }

      override def stop(): Unit = {
        process foreach {
          waitProcess(_, destroy = true)
        }
        process = None
      }

      override def isAlive: Boolean = process.exists { _.isAlive }

      override def logFile: String = cromwellConfiguration.logFile
    }

    JarCromwellProcess(this)
  }
}
