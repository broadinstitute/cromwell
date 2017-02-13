package cromwell.backend.impl.tes

import spray.json._

final case class TesPostResponse(value: String)

final case class TesGetResponse(
                                 jobID: String,
                                 task: TesTaskMessage,
                                 state: String,
                                 logs: Option[Seq[JobLogs]]
                               )


object TesResponseJsonFormatter extends DefaultJsonProtocol {
  implicit val volumeFormat = jsonFormat5(Volume)
  implicit val resourcesFormat = jsonFormat5(Resources)
  implicit val taskParameterFormat = jsonFormat6(TaskParameter)
  implicit val portsFormat = jsonFormat2(Ports)
  implicit val dockerExecutorFormat = jsonFormat7(DockerExecutor)
  implicit val tesTaskMessageFormat = jsonFormat7(TesTaskMessage)
  implicit val jobLogsFormat = jsonFormat8(JobLogs)
  implicit val tesPostResponseFormat = jsonFormat1(TesPostResponse)
  implicit val tesGetResponseFormat = jsonFormat4(TesGetResponse)
}
