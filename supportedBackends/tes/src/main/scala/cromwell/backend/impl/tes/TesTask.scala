package cromwell.backend.impl.tes

import cromwell.backend.BackendJobDescriptor
import cromwell.backend.wdl.OnlyPureFunctions
import wdl4s.values.{WdlFile, WdlSingleFile}

import scala.util.{Failure, Try}


final case class TesTask(jobDescriptor: BackendJobDescriptor) {

  val name = jobDescriptor.key.tag
  // Get the workflow name
  val projectId = jobDescriptor.workflowDescriptor.workflowNamespace.workflow.unqualifiedName
  val description = jobDescriptor.toString
  val taskId = jobDescriptor.toString

  val commands = {
    // The jobDescriptor.inputs are of type WdlValue and need to be casted to WdlFile
    val cliInputs = jobDescriptor.inputs.mapValues(f => WdlFile(f.valueString))

    jobDescriptor
      .call
      .instantiateCommandLine(
        cliInputs,
        OnlyPureFunctions,
        identity
      )
      .map(Seq(_))
  }

  val inputs = jobDescriptor
    .inputs
    .toSeq
    .map {
      case (inputName, f: WdlSingleFile) => TaskParameter(
        Some(inputName),     // Name
        None,           // Description
        Some(f.value),  // Source path
        Some("/tmp"),   // Destination path   NOTE: this must match a Volume.mountPoint in the Resources
        Some("file"),   // Type
        Some(true)      // Create?
      )
    })

//  val outputs = Some(jobDescriptor
//    .call
//    .task
//    .outputs)
  val outputs = None

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

  val dockerImageId = Try(jobDescriptor.runtimeAttributes("docker"))
    .map(_.valueString)
    .recoverWith({
      case e: NoSuchElementException => Failure(ConfigError(
        "Runtime attribute 'docker' is required for TES backend tasks"))
    })

  val dockerExecutor = {
    val workingDirectory = None
    val stdoutPath = Some("/tmp/test_out")
    val stderrPath = None

    for {
      imageId <- dockerImageId
      cmd <- commands
    } yield Seq(DockerExecutor(
      imageId,
      cmd,
      workingDirectory,
      stdoutPath,
      stderrPath))
  }
}


final case class ConfigError(msg: String) extends Exception(msg)

final case class TesTaskMessage(name: String,
                                description: String,
                                projectId: String,
                                taskId: String,
                                inputs: Option[Seq[TaskParameter]],
                                outputs: Option[Seq[TaskParameter]],
                                resources: Resources,
                                docker: Seq[DockerExecutor])

final case class DockerExecutor(imageName: String,
                                cmd: Seq[String],
                                workDir: Option[String],
                                stdout: Option[String],
                                stderr: Option[String])

final case class TaskParameter(name: String,
                               description: Option[String],
                               location: String,
                               path: String,
                               `class`: String,
                               create: Boolean)

final case class Resources(minimumCpuCores: Option[Int],
                           minimumRamGb: Option[Int],
                           preemptible: Option[Boolean],
                           volumes: Option[Seq[Volume]],
                           zones: Option[Seq[String]])

final case class Volume(name: Option[String],
                        sizeGb: Option[Int],
                        source: Option[String],
                        mountPoint: Option[String])
