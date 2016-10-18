package cromwell.backend.impl.tes.util

import spray.json._

final case class TesPostResponse(value: Option[String])

final case class TesGetResponse(jobId: Option[String],
                          task: Option[TesTask],
                          state: Option[String],
                          logs: Option[Seq[Map[String, String]]])

final case class TesTask(name: Option[String],
                   projectId: Option[String],
                   description: Option[String],
                   inputs: Option[Seq[TaskParameter]],
                   outputs: Option[Seq[TaskParameter]],
                   resources: Option[Resources],
                   taskId: Option[String],
                   docker: Option[Seq[DockerExecutor]])

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

object TesResponseJsonFormatter extends DefaultJsonProtocol {
  implicit val volumeFormat = jsonFormat4(Volume)
  implicit val resourcesFormat = jsonFormat5(Resources)
  implicit val taskParameterFormat = jsonFormat6(TaskParameter)
  implicit val dockerExecutorFormat = jsonFormat5(DockerExecutor)
  implicit val tesTaskFormat = jsonFormat8(TesTask)
  implicit val tesPostResponseFormat = jsonFormat1(TesPostResponse)
  implicit val tesGetResponseFormat = jsonFormat4(TesGetResponse)
}
