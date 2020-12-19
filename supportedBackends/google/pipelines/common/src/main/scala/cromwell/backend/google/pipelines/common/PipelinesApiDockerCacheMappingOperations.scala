package cromwell.backend.google.pipelines.common

import _root_.io.circe.parser._
import cats.effect.IO
import com.google.api.services.storage.StorageScopes
import com.google.cloud.storage.{BlobId, Storage, StorageOptions}
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.filesystems.gcs.GcsPathBuilder.ValidFullGcsPath

trait PipelinesApiDockerCacheMappingOperations {

  def generateDockerImageToDiskImageMapping(auth: GoogleAuthMode,
                                            dockerImageCacheManifestFile: ValidFullGcsPath): Map[String, String] = {

    val gcsClient = StorageOptions
      .newBuilder()
      .setCredentials(auth.credentials(Set(StorageScopes.DEVSTORAGE_READ_ONLY)))
      .build
      .getService
    val mappingsFromManifestIO = readDockerImageCacheManifestFileFromGCS(gcsClient, dockerImageCacheManifestFile)
    mappingsFromManifestIO.unsafeRunSync()
  }

  private[common] def readDockerImageCacheManifestFileFromGCS(gcsClient: Storage, gcsPath: ValidFullGcsPath): IO[Map[String, String]] = {
    val manifestFileBlobIo = IO { gcsClient.get(BlobId.of(gcsPath.bucket, gcsPath.path.substring(1))) }
    manifestFileBlobIo flatMap { manifestFileBlob =>
      val jsonStringIo = IO { manifestFileBlob.getContent().map(_.toChar).mkString }
      jsonStringIo.flatMap(jsonStr => IO.fromEither(decode[Map[String, String]](jsonStr)))
    }
  }

}
