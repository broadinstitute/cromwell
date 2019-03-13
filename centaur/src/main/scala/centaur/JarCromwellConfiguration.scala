package centaur

import java.lang.ProcessBuilder.Redirect

import better.files.File
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

      override def displayString: String = command.mkString(" ")

      override def start(): Unit = {
        val processBuilder = new java.lang.ProcessBuilder()
          .command(command: _*)
          .redirectOutput(Redirect.appendTo(File(logFile).toJava))
          .redirectErrorStream(true)
        process = Option(processBuilder.start())
      }

      override def stop(): Unit = {
        process foreach { p =>
          p.getOutputStream.flush()
          p.destroy()
          p.waitFor()
        }
        ()
      }

      override def isAlive: Boolean = process.exists { _.isAlive }

      override def logFile: String = cromwellConfiguration.logFile
    }

    JarCromwellProcess(this)
  }
}
