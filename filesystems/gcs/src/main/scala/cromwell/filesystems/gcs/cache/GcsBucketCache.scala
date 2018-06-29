package cromwell.filesystems.gcs.cache

import java.nio.file.NoSuchFileException

import cats.effect.IO
import com.google.cloud.storage.Storage.{BucketField, BucketGetOption}
import com.google.cloud.storage.{Storage, StorageException}
import com.google.common.cache.Cache
import cromwell.core.path.cache.BucketCache
import cromwell.filesystems.gcs.cache.GcsBucketInformation._

class GcsBucketCache(cloudStorage: Storage, cache: Cache[String, GcsBucketInformation], projectId: String) extends BucketCache[GcsBucketInformation](cache) {
  override protected def retrieve(key: String) = {
    
    def getBucket(bucketGetOptions: BucketGetOption*) = Option(cloudStorage
      .get(key, bucketGetOptions: _*))
      .getOrElse(throw new NoSuchFileException(s"GCS bucket $key does not exist"))

    /*
     * Attempts to get the actual Billing information about the bucket.
     * We need to specify a billing project, because if the bucket does have requester pays, this request
     * will fail without a billing project.
     * Note that even with the project specified, the request can still fail if the credentials used do not have billing
     * permission on that project, in which case we'll still won't know if requester pays is enabled.
     */
    def requesterPaysRequestWithProject: IO[Boolean] = IO {
      val bucket = getBucket(BucketGetOption.userProject(projectId), BucketGetOption.fields(BucketField.BILLING))
      Option(bucket.requesterPays()).exists(_.booleanValue())
    }

    /*
     * Attempts to get metadata about this bucket without specifying a project.
     * If the bucket does have requester pays, the request will fail with a 400 UserProjectMissing:
     * https://cloud.google.com/storage/docs/requester-pays
     */
    def requesterPaysRequestWithoutProject: IO[Boolean] = {
      IO(getBucket(BucketGetOption.fields(BucketField.ID)))
        // If we can get the metadata without specifying a billing project, then requester pays is off
        .map(_ => false)
        .handleErrorWith({
          case storageException: StorageException if isProjectNotProvidedError(storageException) => IO.pure(true)
          case other => IO.raiseError(other)
        })
    }

    // First try with the project
    requesterPaysRequestWithProject
      .handleErrorWith({
        // If it fails because credentials don't have the permission, fallback to trying without project
        case storageException: StorageException if storageException.getCode == 403 &&
          storageException.getMessage.contains(DoesNotHaveServiceUsePermissionError) => requesterPaysRequestWithoutProject
      })
      .map(GcsBucketInformation.apply)
  }
}
