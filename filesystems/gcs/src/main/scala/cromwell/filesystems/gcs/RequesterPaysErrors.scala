package cromwell.filesystems.gcs

import com.google.api.client.googleapis.json.{GoogleJsonError, GoogleJsonResponseException}
import com.google.cloud.storage.StorageException
import org.apache.commons.lang3.StringUtils

object RequesterPaysErrors {
  val BucketIsRequesterPaysErrorCode = 400
  val BucketIsRequesterPaysErrorMessage = "requester pays bucket but no user project"
  val DoesNotHaveServiceUsePermissionErrorCode = 403
  val DoesNotHaveServiceUsePermissionErrorMessage = "does not have serviceusage.services.use"

  def isProjectNotProvidedError(storageException: StorageException) =
    storageException.getCode == BucketIsRequesterPaysErrorCode &&
    StringUtils.contains(storageException.getMessage, BucketIsRequesterPaysErrorMessage)

  def isProjectNotProvidedError(googleJsonError: GoogleJsonError) =
    googleJsonError.getCode == BucketIsRequesterPaysErrorCode &&
      StringUtils.contains(googleJsonError.getMessage, BucketIsRequesterPaysErrorMessage)

  def isProjectNotProvidedError(googleJsonError: GoogleJsonResponseException) =
    googleJsonError.getStatusCode == BucketIsRequesterPaysErrorCode &&
      StringUtils.contains(googleJsonError.getContent, BucketIsRequesterPaysErrorMessage)
}
