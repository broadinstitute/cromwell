package cromwell.backend.impl.tes

import spray.json._

final case class TesPostResponse(value: Option[String])

final case class TesGetResponse(
                                 jobID: Option[String],
                                 task:  Option[TesTaskMessage],
                                 state: Option[String],
                                 logs:  Option[Seq[Map[String, String]]]
                               )


object TesResponseJsonFormatter extends DefaultJsonProtocol {
  implicit val volumeFormat          = jsonFormat4(Volume)
  implicit val resourcesFormat       = jsonFormat5(Resources)
  implicit val taskParameterFormat   = jsonFormat6(TaskParameter)
  implicit val dockerExecutorFormat  = jsonFormat6(DockerExecutor)
  implicit val tesTaskMessageFormat  = jsonFormat7(TesTaskMessage)
  implicit val tesPostResponseFormat = jsonFormat1(TesPostResponse)
  implicit val tesGetResponseFormat  = jsonFormat4(TesGetResponse)
}
