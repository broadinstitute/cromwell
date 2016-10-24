package cromwell.backend.impl.tes.util

import TesTaskCompanion._
import cromwell.backend.BackendJobDescriptor
//import cromwell.backend.wdl.OnlyPureFunctions

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

  def from(jobDescriptor: BackendJobDescriptor): TesTask = {
    // FIXME: this doesn't work - but Seq("echo", "foo") below should be commandLine instead
    //val commandLine = jobDescriptor.call.instantiateCommandLine(Map.empty, OnlyPureFunctions, identity).get

    val docker = jobDescriptor.runtimeAttributes("docker").valueString
    val dockerExecutor = DockerExecutor(Option(docker), Option(Seq("echo", "foo")), None, Option("/tmp/test_out"), None)

    TesTask(Option("TestMD5"),
      Option("My Project"),
      Option("My Desc"),
      None,
      Some(Seq(TaskParameter(None, None, None, Some("/tmp/test_out"), None, None))),
      Some(Resources(None, None, None, Some(Seq(Volume(Some("test_file"), Some(1), None, Some("/tmp")))), None)),
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

