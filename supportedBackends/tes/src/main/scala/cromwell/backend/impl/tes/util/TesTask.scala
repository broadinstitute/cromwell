package cromwell.backend.impl.tes.util

import scala.util.{Failure, Try}
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.wdl.OnlyPureFunctions
import wdl4s.values.{WdlFile, WdlSingleFile}


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

  case class ConfigError(msg: String) extends Exception(msg)

  private def dockerImageId(jobDescriptor: BackendJobDescriptor): Try[String] = {
    Try(jobDescriptor.runtimeAttributes("docker"))
      .map(_.valueString)
      .recoverWith({
        case e: NoSuchElementException => Failure(ConfigError("Runtime attribute 'docker' is required for TES backend tasks"))
      })
  }

  def from(jobDescriptor: BackendJobDescriptor): Try[TesTask] = {

    val taskName = jobDescriptor.key.tag
    val projectId = jobDescriptor.call.parent.get.unqualifiedName
    val description = jobDescriptor.toString
    val taskId = jobDescriptor.toString

    val inputs = Some(jobDescriptor
      .inputs
      .toSeq
      .map {
        case (name, f: WdlSingleFile) => TaskParameter(
          Some(name),     // Name
          Some(name),     // Description
          Some(f.value),  // Source path
          Some("/tmp"),   // Destination path   DEVNOTE: this must match a Volume.mountPoint in the Resources
          Some("file"),   // Type
          Some(true)      // Create?
        )
      })

    val outputs = Some(Seq())

    val volumes = Some(Seq(
      Volume(
        Some("volume name"), // Name
        Some(1),             // Size in GB
        None,                // Source
        Some("/tmp")         // Mount point
      )
    ))

    val resources = Resources(
      None, // Minimum CPU cores
      None, // Minimum RAM in GB
      None, // Preemptible?
      volumes,
      None  // Zones
    )

    for {
      imageId <- dockerImageId(jobDescriptor)
      command <- getCommand(jobDescriptor)
    } yield {
      val dockerExecutor = Seq(DockerExecutor(
        imageId,
        Option(Seq(command)),
        None,
        Option("/tmp/test_out"),
        None))

      TesTask(
        taskName,
        projectId,
        description,
        inputs,
        outputs,
        resources,
        taskId,
        dockerExecutor)
    }
  }

  final case class TesTask(name: String,
                           projectId: String,
                           description: String,
                           inputs: Option[Seq[TaskParameter]],
                           outputs: Option[Seq[TaskParameter]],
                           resources: Resources,
                           taskId: String,
                           docker: Seq[DockerExecutor])

  final case class DockerExecutor(imageName: String,
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
                             minimumRamGb: Option[Int],
                             preemptible: Option[Boolean],
                             volumes: Option[Seq[Volume]],
                             zones: Option[Seq[String]])

  final case class Volume(name: Option[String],
                          sizeGb: Option[Int],
                          source: Option[String],
                          mountPoint: Option[String])
}

