package cromwell.backend.impl.tes.util

import spray.json._

case class TesPostResponse(value: String)

case class TesGetResponse(jobId: String, 
                          task: Map[String, String],
                          state: String, 
                          logs: Map[String, String])

object TesResponseJsonFormatter extends DefaultJsonProtocol {
  implicit val tesPostResponseFormat = jsonFormat1(TesPostResponse)
  implicit val tesGetResponseFormat = jsonFormat4(TesGetResponse)
}
