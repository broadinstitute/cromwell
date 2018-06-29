package cromwell.filesystems.gcs.bucket

import java.io.{IOException, InputStream}
import java.nio.channels.Channels
import java.nio.file.NoSuchFileException

import akka.http.scaladsl.model.{ContentTypes, StatusCodes}
import better.files.File.OpenOptions
import cats.effect.IO
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.cloud.storage.Storage.{BlobSourceOption, BlobTargetOption}
import com.google.cloud.storage.{BlobInfo, Storage, StorageException}
import com.google.common.cache.Cache
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.gcs.cache._
import mouse.all._

import scala.io.Codec

/**
  * Defines behavior on how meta information about buckets can be retrieved and stored.
  */
object GcsBucketInformationPolicies {
  sealed trait GcsBucketInformationPolicy {
    /**
      * Create a reques corresponding to the current policy
      */
    def toRequestHandler(cloudStorage: Storage, projectId: String): GcsRequestHandler
  }

  /**
    * Use to disable access to requester pays buckets
    */
  case object DisabledPolicy extends GcsBucketInformationPolicy {
    override def toRequestHandler(cloudStorage: Storage, projectId: String) = new DisabledRequestHandler(cloudStorage, projectId)
  }

  /**
    * Use to force a request without project set. It is up to the request logic to retry with project set upon failure if appropriate
    */
  case object OnDemandPolicy extends GcsBucketInformationPolicy {
    override def toRequestHandler(cloudStorage: Storage, projectId: String) = new OnDemandRequestHandler(cloudStorage, projectId)
  }

  /**
    * Use to cache bucket information in the provided Guava cache
    */
  case class CachedPolicy(cache: Cache[String, GcsBucketInformation]) extends GcsBucketInformationPolicy {
    override def toRequestHandler(cloudStorage: Storage, projectId: String) = new CachedRequestHandler(cloudStorage, projectId, cache)
  }
}

class DisabledRequestHandler(cloudStorage: Storage, projectId: String) extends GcsRequestHandler(cloudStorage, projectId) {
  override def requesterPays(bucket: String) = RequesterPaysValue.Disabled
}

/**
  * On demand does not lookup nor cache information about buckets. Instead, the request should be tried without project
  * and then retried with project if the error matches.
  */
class OnDemandRequestHandler(cloudStorage: Storage, projectId: String) extends GcsRequestHandler(cloudStorage, projectId) {
  override def requesterPays(bucket: String) = RequesterPaysValue.Unknown
}

/**
  * Cached policy: looks up and caches information about buckets
  */
class CachedRequestHandler(cloudStorage: Storage, projectId: String, cache: Cache[String, GcsBucketInformation]) extends GcsRequestHandler(cloudStorage, projectId) {
  val gcsBucketCache = new GcsBucketCache(cloudStorage, cache, projectId)
  def requesterPays(bucket: String) = RequesterPaysValue.Known(gcsBucketCache.getCachedValue(bucket).unsafeRunSync().requesterPays)
}

/**
  * Implements read and write requests according to the GcsBucketInformationPolicy
  */
abstract class GcsRequestHandler(cloudStorage: Storage, projectId: String) {
  
  /**
    * Write content to the provided path.
    * The request might be executed a second time with a billing project set if the bucket has RP and depending on the policy.
    */
  def write(path: GcsPath, content: String, openOptions: OpenOptions, codec: Codec): IO[GcsPath] = recoverFromProjectNotProvided(path, "write to") {
    _write(path, content, openOptions, codec)
  }

  /**
    * Return an input stream for the provided path.
    * The request might be executed a second time with a billing project set if the bucket has RP and depending on the policy.
    */
  def inputStream(path: GcsPath): IO[InputStream] = recoverFromProjectNotProvided(path, "read from") {
    _inputStream(path)
  }
  
  private def _write(path: GcsPath, content: String, openOptions: OpenOptions, codec: Codec)(requesterPaysValue: RequesterPaysValue): GcsPath = {
    cloudStorage.create(
      BlobInfo.newBuilder(path.blob)
        .setContentType(ContentTypes.`text/plain(UTF-8)`.value)
        .build(),
      content.getBytes(codec.charSet),
      requesterPaysValue.withProject.option(userProjectBlobTarget).toList.flatten: _*
    )
    path
  }

  private def _inputStream(path: GcsPath)(requesterPaysValue: RequesterPaysValue): InputStream = {
    Channels.newInputStream(cloudStorage.reader(path.blob, requesterPaysValue.withProject.option(userProjectBlobSource).toList.flatten: _*))
  }

  def requesterPays(bucket: String): RequesterPaysValue

  // If the request fails because no project was passed, recover the request, this time setting the project
  private def recoverFromProjectNotProvided[A](path: GcsPath, action: String)(f: RequesterPaysValue => A) = {
    val requesterPaysValue = requesterPays(path.blob.getBucket)

    IO(f(requesterPaysValue))
      .handleErrorWith({
        // Only retry with the the project if the error is right and requester pays is not disabled
        case error: StorageException if GcsBucketInformation.isProjectNotProvidedError(error) && requesterPaysValue.enabled => 
          IO(f(RequesterPaysValue.Known(true)))
        // Wrap file not founds in NoSuchFileException for better error reporting
        case e: GoogleJsonResponseException if e.getStatusCode == StatusCodes.NotFound.intValue =>
          IO.raiseError(new NoSuchFileException(path.pathAsString))
        // Wrap error in an IOException containing the name of the file  
        case e => 
          IO.raiseError(new IOException(s"Could not access $action ${path.pathAsString}: ${e.getMessage}", e))
      })
  }

  private val userProjectBlobTarget: List[BlobTargetOption] = List(BlobTargetOption.userProject(projectId))
  private val userProjectBlobSource: List[BlobSourceOption] = List(BlobSourceOption.userProject(projectId))
}
