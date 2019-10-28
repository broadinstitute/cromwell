package cromwell.services

import cromwell.services.metadata.MetadataService.{BuildMetadataJsonAction, MetadataServiceResponse}
import spray.json.JsObject

package object metadata {
  type QueryParameters = Seq[QueryParameter]
}

sealed trait MetadataJsonResponse extends MetadataServiceResponse { def originalRequest: BuildMetadataJsonAction }
final case class SuccessfulMetadataJsonResponse(originalRequest: BuildMetadataJsonAction, responseJson: JsObject) extends MetadataJsonResponse
final case class FailedMetadataJsonResponse(originalRequest: BuildMetadataJsonAction, reason: Throwable) extends MetadataJsonResponse
