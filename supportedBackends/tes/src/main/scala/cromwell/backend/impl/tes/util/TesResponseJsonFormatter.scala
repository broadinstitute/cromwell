package cromwell.backend.impl.tes.util

import cromwell.backend.impl.tes.util.TesTaskCompanion.{DockerExecutor, Resources, TaskParameter, Volume, TesTask}
import spray.json._

final case class TesPostResponse(value: Option[String])

final case class TesGetResponse(jobId: Option[String],
                          task: Option[TesTask],
                          state: Option[String],
                          logs: Option[Seq[Map[String, String]]])


object TesResponseJsonFormatter extends DefaultJsonProtocol {
  implicit val volumeFormat = jsonFormat4(Volume)
  implicit val resourcesFormat = jsonFormat5(Resources)
  implicit val taskParameterFormat = jsonFormat6(TaskParameter)
  implicit val dockerExecutorFormat = jsonFormat5(DockerExecutor)
  implicit val tesTaskFormat = jsonFormat8(TesTask)
  implicit val tesPostResponseFormat = jsonFormat1(TesPostResponse)
  implicit val tesGetResponseFormat = jsonFormat4(TesGetResponse)
}