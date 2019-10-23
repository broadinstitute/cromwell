package cromwell.services

import cromwell.services.metadata.MetadataService.{MetadataReadAction, MetadataServiceResponse}
import spray.json.JsObject

package object metadata {
  type QueryParameters = Seq[QueryParameter]
}

sealed trait BuildMetadataResponse extends MetadataServiceResponse { def originalRequest: MetadataReadAction }
final case class BuiltMetadataResponse(originalRequest: MetadataReadAction, responseJson: JsObject) extends BuildMetadataResponse
final case class FailedMetadataResponse(originalRequest: MetadataReadAction, reason: Throwable) extends BuildMetadataResponse
