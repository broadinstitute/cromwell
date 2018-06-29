package cromwell.filesystems.gcs.cache

import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.cloud.storage.StorageException

object GcsBucketInformation {
  val BucketIsRequesterPaysErrorCode = 400
  val BucketIsRequesterPaysErrorMessage = "Bucket is requester pays bucket but no user project provided."
  val DoesNotHaveServiceUsePermissionError = "does not have serviceusage.services.use"

  def isProjectNotProvidedError(storageException: StorageException) = 
    storageException.getCode == BucketIsRequesterPaysErrorCode &&
    storageException.getMessage == BucketIsRequesterPaysErrorMessage

  def isProjectNotProvidedError(googleJsonError: GoogleJsonError) =
    googleJsonError.getCode == BucketIsRequesterPaysErrorCode &&
      googleJsonError.getMessage == BucketIsRequesterPaysErrorMessage
}

case class GcsBucketInformation(requesterPays: Boolean)
