package cromwell.filesystems.gcs.cache

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
import cromwell.filesystems.gcs.cache.GcsRequestHandler.RequesterPaysValue
import mouse.all._

import scala.io.Codec
import scala.util.{Failure, Success, Try}

object GcsRequestHandler {
  sealed trait RequesterPaysValue {
    def withProject: Boolean
    def enabled: Boolean
  }
  object RequesterPaysValue{
    case class Known(value: Boolean) extends RequesterPaysValue {
      override def withProject = value
      override def enabled = true
    }
    case object Unknown extends RequesterPaysValue {
      override def withProject = false
      override def enabled = true
    }
    case object Disabled extends RequesterPaysValue {
      override def withProject = false
      override def enabled = false
    }
  }
}

/**
  * Implements read and write requests according to the GcsBucketInformationPolicy
  */
abstract class GcsRequestHandler(cloudStorage: Storage, projectId: String) {
  def write(path: GcsPath, content: String, openOptions: OpenOptions, codec: Codec): IO[GcsPath] = recoverFromProjectNotProvided(path) {
    _write(path, content, openOptions, codec)
  }

  def inputStream(path: GcsPath): IO[InputStream] = recoverFromProjectNotProvided(path) {
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
    Try{
      Channels.newInputStream(cloudStorage.reader(path.blob, requesterPaysValue.withProject.option(userProjectBlobSource).toList.flatten: _*))
    } match {
      case Success(inputStream) => inputStream
      case Failure(e: GoogleJsonResponseException) if e.getStatusCode == StatusCodes.NotFound.intValue =>
        throw new NoSuchFileException(path.pathAsString)
      case Failure(e) => e.getMessage
        throw new IOException(s"Failed to open an input stream for ${path.pathAsString}: ${e.getMessage}", e)
    }
  }

  def requesterPays(bucket: String): RequesterPaysValue

  // If the request fails because no project was passed, recover the request, this time setting the project
  private def recoverFromProjectNotProvided[A](path: GcsPath)(f: RequesterPaysValue => A) = {
    val requesterPaysValue = requesterPays(path.blob.getBucket)
    IO(f(requesterPaysValue))
      .handleErrorWith({
        // Only retry with the the project if the error is right and requester pays is not disabled
        case error: StorageException if GcsBucketInformation.isProjectNotProvidedError(error) &&
          requesterPaysValue.enabled => IO(f(RequesterPaysValue.Known(true)))
      })
  }

  private val userProjectBlobTarget: List[BlobTargetOption] = List(BlobTargetOption.userProject(projectId))
  private val userProjectBlobSource: List[BlobSourceOption] = List(BlobSourceOption.userProject(projectId))
}

class DisabledRequestHandler(cloudStorage: Storage, projectId: String) extends GcsRequestHandler(cloudStorage, projectId) {
  override def requesterPays(bucket: String) = RequesterPaysValue.Disabled
}

/**
  * On demand does not lookup nor cache information about buckets.
  * Therefore it always returns None as it never knows if a requester has requester pays or not
  */
class OnDemandRequestHandler(cloudStorage: Storage, projectId: String) extends GcsRequestHandler(cloudStorage, projectId) {
  override def requesterPays(bucket: String) = RequesterPaysValue.Unknown
}

/**
  * Cached policy looks up and caches information about buckets
  */
class CachedRequestHandler(cloudStorage: Storage, projectId: String, cache: Cache[String, GcsBucketInformation]) extends GcsRequestHandler(cloudStorage, projectId) {
  val gcsBucketCache = new GcsBucketCache(cloudStorage, cache, projectId)
  def requesterPays(bucket: String) = RequesterPaysValue.Known(gcsBucketCache.getCachedValue(bucket).unsafeRunSync().requesterPays)
}

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
 
