package cromwell.services

import cromwell.core.WorkflowId
import cromwell.services.metadata.MetadataService.{BuildMetadataJsonAction, MetadataServiceResponse}
import spray.json.JsObject

package object metadata {
  type QueryParameters = Seq[QueryParameter]
}

sealed trait MetadataJsonResponse extends MetadataServiceResponse { def originalRequest: BuildMetadataJsonAction }
final case class SuccessfulMetadataJsonResponse(originalRequest: BuildMetadataJsonAction, responseJson: JsObject) extends MetadataJsonResponse
final case class FailedMetadataJsonResponse(originalRequest: BuildMetadataJsonAction, reason: Throwable) extends MetadataJsonResponse

final class MetadataTooLargeException(workflowId: WorkflowId,
                                metadataSizeRows: Option[Int],
                                metadataLimitRows: Int)
  extends RuntimeException(s"Metadata for workflow $workflowId exists in" +
    s"database, but cannot be served. This is done in order to avoid Cromwell failure: metadata is too large - " +
    s"$metadataSizeRows rows, and may cause Cromwell instance to die on attempt to read it in memory. Configured " +
    s"metadata safety limit is $metadataLimitRows")
