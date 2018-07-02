package cromwell.filesystems.gcs

import com.google.api.client.googleapis.json.{GoogleJsonError, GoogleJsonResponseException}
import com.google.cloud.storage.StorageException

object RequesterPaysErrors {
  val BucketIsRequesterPaysErrorCode = 400
  val BucketIsRequesterPaysErrorMessage = "Bucket is requester pays bucket but no user project provided."
  val DoesNotHaveServiceUsePermissionErrorCode = 403
  val DoesNotHaveServiceUsePermissionErrorMessage = "does not have serviceusage.services.use"

  def isProjectNotProvidedError(storageException: StorageException) = 
    storageException.getCode == BucketIsRequesterPaysErrorCode &&
    storageException.getMessage == BucketIsRequesterPaysErrorMessage

  def isProjectNotProvidedError(googleJsonError: GoogleJsonError) =
    googleJsonError.getCode == BucketIsRequesterPaysErrorCode &&
      googleJsonError.getMessage == BucketIsRequesterPaysErrorMessage

  def isProjectNotProvidedError(googleJsonError: GoogleJsonResponseException) =
    googleJsonError.getStatusCode == BucketIsRequesterPaysErrorCode &&
      googleJsonError.getContent == BucketIsRequesterPaysErrorMessage
}
