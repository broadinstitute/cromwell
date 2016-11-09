package cromwell.backend.impl.tes.util

import TesTaskCompanion._
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.wdl.OnlyPureFunctions
import wdl4s.values.{WdlFile, WdlSingleFile}

import scala.util.Try

final case class TesTask(name: Option[String],
                         projectId: Option[String],
                         description: Option[String],
                         inputs: Option[Seq[TaskParameter]],
                         outputs: Option[Seq[TaskParameter]],
                         resources: Option[Resources],
                         taskId: Option[String],
                         docker: Option[Seq[DockerExecutor]])

/*
  Because of spray serialization support we can't have the companion object named TesTask, and that also has
  an impact on what should be `apply` not being appropriate as well
 */
object TesTaskCompanion {

  private def getCommand(jobDescriptor: BackendJobDescriptor): Try[String] = {
    // The jobDescriptor.inputs are of type WdlValue and need to be casted to WdlFile
    val cliInputs = jobDescriptor.inputs.mapValues(f => WdlFile(f.valueString))

    jobDescriptor
      .call
      .instantiateCommandLine(
        cliInputs,
        OnlyPureFunctions,
        identity
      )
  }

  def from(jobDescriptor: BackendJobDescriptor): TesTask = {

    // Get the "command" section of the WDL task.
    val command = getCommand(jobDescriptor).get

    // TODO it's possible that the "docker" key doesn't exist
    val docker = jobDescriptor
      .runtimeAttributes("docker")
      .valueString

    val dockerExecutor = DockerExecutor(
      Option(docker),
      Option(Seq(command)),
      None,
      Option("/tmp/test_out"),
      None)

    val inputs = jobDescriptor
      .inputs
      .toSeq
      .map {
        case (lqName, f: WdlSingleFile) => TaskParameter(
          Some(lqName),
          Some("description"),
          Some(f.value),
          Some("path"),
          Some("file"),
          Some(true)
        )
      }

    val resources = Resources(
      None,
      None,
      None,
      Some(
        Seq(
          Volume(
            Some("test_file"),
            Some(1),
            None,
            Some("/tmp")
          ))),
      None
    )

    TesTask(
      Option("TestMD5"),
      Option("My Project"),
      Option("My Desc"),
      Some(inputs),
      None,
      Some(resources),
      None,
      Option(Seq(dockerExecutor)))
  }

  final case class DockerExecutor(imageName: Option[String],
                                  cmd: Option[Seq[String]],
                                  workDir: Option[String],
                                  stdout: Option[String],
                                  stderr: Option[String])

  final case class TaskParameter(name: Option[String],
                                 description: Option[String],
                                 location: Option[String],
                                 path: Option[String],
                                 `class`: Option[String],
                                 create: Option[Boolean])

  final case class Resources(minimumCpuCores: Option[Int],
                             preemtible: Option[Boolean],
                             minimumRamGb: Option[Int],
                             volumes: Option[Seq[Volume]],
                             zones: Option[Seq[String]])

  final case class Volume(name: Option[String],
                          sizeGb: Option[Int],
                          source: Option[String],
                          mountPoint: Option[String])
}

