package cromwell.backend.impl.sfs.config

import com.typesafe.config.Config
import cromwell.backend.impl.sfs.config.ConfigConstants._
import net.ceedubs.ficus.Ficus._
import wdl.draft2.model.{WdlNamespace, WdlTask}
import wom.core.WorkflowSource

import scala.util.{Failure, Success}

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

  private val killDockerCommandOption = backendConfig.as[Option[String]](KillDockerConfig)
  private val killDockerSourceOption = killDockerCommandOption.map(makeWdlSource(KillDockerTask, _, jobIdRuntimeAttributes))

  private val checkAliveCommandOption = backendConfig.as[Option[String]](CheckAliveConfig)
  private val checkAliveSourceOption = checkAliveCommandOption.map(makeWdlSource(
    CheckAliveTask, _, jobIdRuntimeAttributes))

  private val workflowSource =
    s"""
       |${submitSourceOption getOrElse ""}
       |${submitDockerSourceOption getOrElse ""}
       |${killSourceOption getOrElse ""}
       |${killDockerSourceOption getOrElse ""}
       |${checkAliveSourceOption getOrElse ""}
       |""".stripMargin.trim

  /**
    * The wdl namespace containing the submit, kill, and check alive tasks.
    */
  val wdlNamespace = {
    WdlNamespace.loadUsingSource(workflowSource, None, None) match {
      case Success(ns) => ns
      case Failure(f) => throw new RuntimeException(s"Error parsing generated wdl:\n$workflowSource".stripMargin, f)
    }
  }

  private val runtimeAttributesTask = makeTask(RuntimeAttributesTask, "", configRuntimeAttributes)

  /**
    * The declarations of runtime attributes.
    */
  val runtimeDeclarations = runtimeAttributesTask.declarations

}

object ConfigWdlNamespace {
  private def makeWdlSource(taskName: String, command: String, declarations: String): WorkflowSource = {
    s"""
       |task $taskName {
       |$declarations
       |command {
       |$command
       |}
       |}
       |""".stripMargin
  }

  private def makeTask(taskName: String, command: String, declarations: String): WdlTask = {
    val workflowSource = makeWdlSource(taskName, command, declarations)
    val wdlNamespace = WdlNamespace.loadUsingSource(workflowSource, None, None).get
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
       |String $JobShellInput
       |""".stripMargin

  /**
    * Extra inputs only filled in for submit-docker.
    */
  private val submitDockerRuntimeAttributes =
    s"""
       |String $DockerCwdInput
       |String $DockerCidInput
       |String $DockerScriptInput
       |String $DockerStdoutInput
       |String $DockerStderrInput
       |""".stripMargin

  /**
    * Inputs for the kill and check alive commands.
    */
  private val jobIdRuntimeAttributes =
    s"""
       |String $JobIdInput
       |String $DockerCidInput
       |String $JobShellInput
       |""".stripMargin
}
