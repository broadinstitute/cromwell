package cromwell.api.model

import spray.json.DefaultJsonProtocol

object CromwellFailedResponseExceptionJsonSupport extends DefaultJsonProtocol {
  implicit val CromwellFailedResponseExceptionFormat = jsonFormat2(CromwellFailedResponseException)
}

case class CromwellFailedResponseException(status: String, message: String) extends Exception(message)
