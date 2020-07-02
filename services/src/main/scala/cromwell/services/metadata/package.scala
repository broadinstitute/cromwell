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

class MetadataTooLargeException(message: String) extends RuntimeException(message) {
  // we don't want to show stacktrace for these errors, so suppress stacktrace's population
  override def fillInStackTrace(): Throwable = this
}

final class MetadataTooLargeNumberOfRowsException(workflowId: WorkflowId, metadataSizeRows: Int, metadataLimitRows: Int)
  extends MetadataTooLargeException(s"Metadata for workflow $workflowId exists in " +
    s"database, but cannot be served. This is done in order to avoid Cromwell failure: metadata is too large - " +
    s"$metadataSizeRows rows, and may cause Cromwell instance to die on attempt to read it in memory. Configured " +
    s"metadata safety limit is $metadataLimitRows.")

final class MetadataTooLargeTimeoutException(workflowId: WorkflowId)
  extends MetadataTooLargeException(s"Metadata for workflow $workflowId exists in " +
    s"database, but cannot be served. This is done in order to avoid Cromwell failure: metadata is probably too " +
    s"large - timeout occurred on attempt to fetch it from the database.")
