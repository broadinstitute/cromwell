package cromwell.filesystems.gcs

import com.google.api.client.googleapis.json.{GoogleJsonError, GoogleJsonResponseException}
import com.google.cloud.storage.StorageException

object RequesterPaysErrors {
  val BucketIsRequesterPaysErrorCode = 400
  val BucketIsRequesterPaysErrorMessage = "Bucket is requester pays bucket but no user project provided."
  val DoesNotHaveServiceUsePermissionErrorCode = 403
  val DoesNotHaveServiceUsePermissionErrorMessage = "does not have serviceusage.services.use"

  def retryWithBillingProject(storageException: StorageException) = isProjectNotProvidedError || isServiceUsageNotProvided

  def retryWithBillingProject(googleJsonError: GoogleJsonError) = isProjectNotProvidedError || isServiceUsageNotProvided

  def retryWithBillingProject(googleJsonError: GoogleJsonResponseException) = isProjectNotProvidedError || isServiceUsageNotProvided

  def isProjectNotProvidedError(storageException: StorageException) = 
    storageException.getCode == BucketIsRequesterPaysErrorCode &&
    storageException.getMessage == BucketIsRequesterPaysErrorMessage

  def isProjectNotProvidedError(googleJsonError: GoogleJsonError) =
    googleJsonError.getCode == BucketIsRequesterPaysErrorCode &&
      googleJsonError.getMessage == BucketIsRequesterPaysErrorMessage

  def isProjectNotProvidedError(googleJsonError: GoogleJsonResponseException) =
    googleJsonError.getStatusCode == BucketIsRequesterPaysErrorCode &&
      googleJsonError.getContent == BucketIsRequesterPaysErrorMessage

  def isServiceUsageNotProvided(storageException: StorageException) =
    storageException.getCode == DoesNotHaveServiceUsePermissionErrorCode &&
      storageException.getMessage == DoesNotHaveServiceUsePermissionErrorMessage

  def isServiceUsageNotProvided(googleJsonError: GoogleJsonError) =
    googleJsonError.getCode == DoesNotHaveServiceUsePermissionErrorCode &&
      googleJsonError.getMessage == DoesNotHaveServiceUsePermissionErrorMessage

  def isServiceUsageNotProvided(googleJsonError: GoogleJsonResponseException) =
    googleJsonError.getStatusCode == DoesNotHaveServiceUsePermissionErrorCode &&
      googleJsonError.getContent == DoesNotHaveServiceUsePermissionErrorMessage
}
