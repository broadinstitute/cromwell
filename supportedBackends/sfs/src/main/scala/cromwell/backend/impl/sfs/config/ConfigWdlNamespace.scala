package cromwell.backend.impl.sfs.config

import com.typesafe.config.Config
import cromwell.backend.impl.sfs.config.ConfigConstants._
import net.ceedubs.ficus.Ficus._
import wdl4s._

/**
  * Builds a wdl namespace from the config.
  *
  * @param backendConfig The config section for a backend instance.
  */
class ConfigWdlNamespace(backendConfig: Config) {

  import ConfigWdlNamespace._

  private val configRuntimeAttributes = backendConfig.as[Option[String]](RuntimeAttributesConfig).getOrElse("")

  private val submitCommandOption = backendConfig.as[Option[String]](SubmitConfig)
  private val submitSourceOption = submitCommandOption.map(makeWdlSource(
    SubmitTask, _, submitRuntimeAttributes + configRuntimeAttributes))

  private val submitDockerCommandOption = backendConfig.as[Option[String]](SubmitDockerConfig)
  private val submitDockerSourceOption = submitDockerCommandOption.map(makeWdlSource(
    SubmitDockerTask, _, submitRuntimeAttributes + submitDockerRuntimeAttributes + configRuntimeAttributes))

  private val killCommandOption = backendConfig.as[Option[String]](KillConfig)
  private val killSourceOption = killCommandOption.map(makeWdlSource(KillTask, _, jobIdRuntimeAttributes))

  private val checkAliveCommandOption = backendConfig.as[Option[String]](CheckAliveConfig)
  private val checkAliveSourceOption = checkAliveCommandOption.map(makeWdlSource(
    CheckAliveTask, _, jobIdRuntimeAttributes))

  private val wdlSource =
    s"""
       |${submitSourceOption getOrElse ""}
       |${submitDockerSourceOption getOrElse ""}
       |${killSourceOption getOrElse ""}
       |${checkAliveSourceOption getOrElse ""}
       |""".stripMargin.trim

  /**
    * The wdl namespace containing the submit, kill, and check alive tasks.
    */
  val wdlNamespace = {
    try {
      WdlNamespace.loadUsingSource(wdlSource, None, None)
    } catch {
      case exception: Exception =>
        throw new RuntimeException(s"Error parsing generated wdl:\n$wdlSource".stripMargin, exception)
    }
  }

  private val runtimeAttributesTask = makeTask(RuntimeAttributesTask, "", configRuntimeAttributes)

  /**
    * The declarations of runtime attributes.
    */
  val runtimeDeclarations = runtimeAttributesTask.declarations

}

object ConfigWdlNamespace {
  private def makeWdlSource(taskName: String, command: String, declarations: String): WdlSource = {
    s"""
       |task $taskName {
       |$declarations
       |command {
       |$command
       |}
       |}
       |""".stripMargin
  }

  private def makeTask(taskName: String, command: String, declarations: String): Task = {
    val wdlSource = makeWdlSource(taskName, command, declarations)
    val wdlNamespace = WdlNamespace.loadUsingSource(wdlSource, None, None)
    wdlNamespace.findTask(taskName).getOrElse(throw new RuntimeException(s"Couldn't find task $taskName"))
  }

  /**
    * Extra inputs that will be filled in for both submit and submit-docker.
    */
  private val submitRuntimeAttributes =
    s"""
       |String $JobIdInput
       |String $JobNameInput
       |String $CwdInput
       |String $StdoutInput
       |String $StderrInput
       |String $ScriptInput
       |""".stripMargin

  /**
    * Extra inputs only filled in for submit-docker.
    */
  private val submitDockerRuntimeAttributes =
    s"""
       |String $DockerCwdInput
       |""".stripMargin

  /**
    * Inputs for the kill and check alive commands.
    */
  private val jobIdRuntimeAttributes =
    s"""
       |String $JobIdInput
       |""".stripMargin
}
