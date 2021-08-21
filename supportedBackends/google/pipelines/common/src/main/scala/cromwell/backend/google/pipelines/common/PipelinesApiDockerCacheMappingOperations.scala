package cromwell.backend.google.pipelines.common

import _root_.io.circe.generic.auto._
import _root_.io.circe.parser._
import cats.effect.IO
import com.google.api.services.storage.StorageScopes
import com.google.cloud.storage.{BlobId, Storage, StorageOptions}
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.logging.JobLogger
import cromwell.filesystems.gcs.GcsPathBuilder.ValidFullGcsPath

import scala.util.control.NoStackTrace

case class DockerImageCacheEntry(dockerImageDigest: String, diskImageName: String)
case class DockerImageCacheManifest(manifestFormatVersion: Int, dockerImageCacheMap: Map[String, DockerImageCacheEntry])

trait PipelinesApiDockerCacheMappingOperations {

  private val CURRENT_SUPPORTED_MANIFEST_FORMAT_VERSION = 2
  private class DockerImageManifestVersionError(message: String) extends RuntimeException(message) with NoStackTrace

  def generateDockerImageToDiskImageMapping(auth: GoogleAuthMode,
                                            dockerImageCacheManifestFile: ValidFullGcsPath): Map[String, DockerImageCacheEntry] = {

    val gcsClient = StorageOptions
      .newBuilder()
      .setCredentials(auth.credentials(Set(StorageScopes.DEVSTORAGE_READ_ONLY)))
      .build
      .getService
    val mappingsFromManifestIO = readDockerImageCacheManifestFileFromGCS(gcsClient, dockerImageCacheManifestFile)
    mappingsFromManifestIO.map(_.dockerImageCacheMap).unsafeRunSync()
  }

  def getDockerCacheDiskImageForAJob(dockerImageToCacheDiskImageMappingOpt: Option[Map[String, DockerImageCacheEntry]],
                                     dockerImageAsSpecifiedByUser: String,
                                     dockerImageWithDigest: String,
                                     jobLogger: JobLogger): Option[String] = {
    dockerImageToCacheDiskImageMappingOpt
      .flatMap(_.get(dockerImageAsSpecifiedByUser))
      .filter { cachedDockerImageDigestAndDiskName =>
        val hashStartingPositionInActualDockerImage = dockerImageWithDigest.indexOf('@')
        if (hashStartingPositionInActualDockerImage != -1) {
          val actualDigestOfDesiredDockerImage = dockerImageWithDigest.substring(hashStartingPositionInActualDockerImage + 1)
          if (cachedDockerImageDigestAndDiskName.dockerImageDigest == actualDigestOfDesiredDockerImage) {
            true
          } else {
            jobLogger.info(s"Cached Docker image digest mismatch. Requested docker image $dockerImageAsSpecifiedByUser has different digest than " +
              s"corresponding cached image located at the ${cachedDockerImageDigestAndDiskName.diskImageName} disk image. " +
              s"Digest of requested image is $actualDigestOfDesiredDockerImage, but digest of cached image is ${cachedDockerImageDigestAndDiskName.dockerImageDigest}. " +
              s"Docker image cache feature will not be used for this task.")

            false
          }
        } else {
          jobLogger.error(s"Programmer error ! Odd docker image name where supposed to be name with digest: $dockerImageWithDigest")
          false
        }
      }
      .map(_.diskImageName)
  }

  private[common] def readDockerImageCacheManifestFileFromGCS(gcsClient: Storage, gcsPath: ValidFullGcsPath): IO[DockerImageCacheManifest] = {
    val manifestFileBlobIo = IO { gcsClient.get(BlobId.of(gcsPath.bucket, gcsPath.path.substring(1))) }
    manifestFileBlobIo flatMap { manifestFileBlob =>
      val jsonStringIo = IO { manifestFileBlob.getContent().map(_.toChar).mkString }
      jsonStringIo.flatMap { jsonStr =>
        decode[DockerImageCacheManifest](jsonStr) match {
          case Left(error) => IO.raiseError(error)
          case Right(parsedManifest) =>
            if (parsedManifest.manifestFormatVersion == CURRENT_SUPPORTED_MANIFEST_FORMAT_VERSION) {
              IO.pure(parsedManifest)
            } else {
              IO.raiseError(new DockerImageManifestVersionError(s"Current supported docker image cache manifest format version " +
                s"is $CURRENT_SUPPORTED_MANIFEST_FORMAT_VERSION, but got ${parsedManifest.manifestFormatVersion}"))
            }
        }
      }
    }
  }

}
