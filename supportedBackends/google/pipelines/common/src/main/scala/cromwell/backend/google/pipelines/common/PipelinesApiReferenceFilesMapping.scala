package cromwell.backend.google.pipelines.common

import java.util

import cats.effect.IO
import cats.implicits._
import com.google.api.services.storage.StorageScopes
import com.google.cloud.storage.{BlobId, Storage, StorageOptions}
import com.google.common.io.BaseEncoding
import com.google.common.primitives.Longs
import cromwell.backend.google.pipelines.common.io.PipelinesApiReferenceFilesDisk
import cromwell.filesystems.gcs.GcsPathBuilder
import cromwell.filesystems.gcs.GcsPathBuilder.ValidFullGcsPath
import PipelinesApiReferenceFilesMapping._
import cromwell.cloudsupport.gcp.auth.{ApplicationDefaultMode, GoogleAuthMode}
import org.slf4j.{Logger, LoggerFactory}

/**
 * This class contains logic related to public Broad reference files stored in public GCS bucket.
 * During instance creation it will read and parse Broad reference disk manifest files from GCS and populate internal
 * map containing relations between reference file paths and reference disks. This may take significant time, so the
 * best way to do it is during Cromwell startup.
 *
 * @param auth - authentication to be use to access GCS
 * @param referenceDiskLocalizationManifestFiles - optional list of reference disk manifest files to parse
 */
class PipelinesApiReferenceFilesMapping(auth: GoogleAuthMode,
                                        referenceDiskLocalizationManifestFiles: Option[List[ValidFullGcsPath]]) {

  private val logger: Logger = LoggerFactory.getLogger(getClass)

  private val validReferenceFilesMap: Map[String, PipelinesApiReferenceFilesDisk] = {
    val validReferenceFilesMapIO = referenceDiskLocalizationManifestFiles match {
      case None =>
        IO.pure(Map.empty[String, PipelinesApiReferenceFilesDisk])
      case Some(manifestFilesPaths) =>
        val manifestFilesIo = manifestFilesPaths traverse { manifestFilePath => readReferenceDiskManifestFileFromGCS(manifestFilePath) }
        manifestFilesIo flatMap { manifestFiles =>
          manifestFiles.traverse(manifestFileToMapOfValidReferences).map(_.flatten.toMap)
        }
    }
    validReferenceFilesMapIO.unsafeRunSync()
  }

  def getReferenceDisksToMount(inputFilePaths: Set[String]): List[PipelinesApiReferenceFilesDisk] = {
    validReferenceFilesMap.filterKeys(key => inputFilePaths.contains(s"gs://$key")).values.toList.distinct
  }

  private def createGcsClient(): Storage = StorageOptions
    .newBuilder()
    .setCredentials(auth.credentials(Set(StorageScopes.DEVSTORAGE_READ_ONLY)))
    .build
    .getService

  protected def readReferenceDiskManifestFileFromGCS(gcsPath: ValidFullGcsPath): IO[ManifestFile] = {
    import spray.json._
    object ManifestFileJsonProtocol extends DefaultJsonProtocol {
      implicit val referenceFileFormat: RootJsonFormat[ReferenceFile] = jsonFormat2(ReferenceFile)
      implicit val manifestFileFormat: RootJsonFormat[ManifestFile] = jsonFormat3(ManifestFile)
    }
    import ManifestFileJsonProtocol._

    val gcsClient = createGcsClient()
    val manifestFileBlobIo = IO { gcsClient.get(BlobId.of(gcsPath.bucket, gcsPath.path.substring(1))) }
    manifestFileBlobIo flatMap { manifestFileBlob =>
      val jsonStringIo = IO { manifestFileBlob.getContent().map(_.toChar).mkString }
      jsonStringIo.map(_.parseJson.convertTo[ManifestFile])
    }
  }

  protected def isReferenceFileChecksumValid(gcsClient: Storage, referenceFile: ReferenceFile): IO[Boolean] = {
    GcsPathBuilder.validateGcsPath(s"gs://${referenceFile.path}") match {
      case ValidFullGcsPath(bucket, path) =>
        val gcsFileBlobIo = IO { gcsClient.get(BlobId.of(bucket, path.substring(1))) }
        val crc32cFromManifest = BaseEncoding.base64.encode(
          // drop 4 leading bytes from Long crc32c value
          // https://stackoverflow.com/a/25111119/1794750
          util.Arrays.copyOfRange(Longs.toByteArray(referenceFile.crc32c), 4, 8)
        )
        gcsFileBlobIo.map(gcsFileBlob => gcsFileBlob.getCrc32c === crc32cFromManifest)
      case _ =>
        IO.pure(false)
    }
  }

  private def manifestFileToMapOfValidReferences(manifestFile: ManifestFile): IO[Map[String, PipelinesApiReferenceFilesDisk]] = {
    val gcsClient = createGcsClient()
    val refDisk = PipelinesApiReferenceFilesDisk(manifestFile.imageIdentifier, manifestFile.diskSizeGb)
    val allReferenceFilesFromManifestMap = manifestFile.files.map(refFile => (refFile, refDisk)).toMap

    // validate checksums of reference files and keep only matching ones
    val validReferenceFilesFromManifestMapFuture = for {
      checkedFiles <-
        allReferenceFilesFromManifestMap
          .keySet
          .map(refFile =>
            isReferenceFileChecksumValid(gcsClient, refFile).map(isValid => (refFile, isValid))
          ).toList.sequence
    } yield allReferenceFilesFromManifestMap.filterKeys { key =>
      checkedFiles.collect {
        case (refFile, isValid) if isValid => refFile
      }.contains(key)
    }

    validReferenceFilesFromManifestMapFuture map { validReferenceFilesFromManifestMap =>
      val invalidReferenceFiles = allReferenceFilesFromManifestMap.keySet -- validReferenceFilesFromManifestMap.keySet
      if (invalidReferenceFiles.nonEmpty) {
        logger.warn(s"The following files listed in references manifest have checksum mismatch with actual files in GCS: ${invalidReferenceFiles.mkString(",")}")
      }
      validReferenceFilesFromManifestMap map {
        case (refFile, disk) => (refFile.path, disk)
      }
    }
  }
}

object PipelinesApiReferenceFilesMapping {
  case class ReferenceFile(path: String, crc32c: Long)
  case class ManifestFile(imageIdentifier: String, diskSizeGb: Int, files: List[ReferenceFile])

  def apply(auth: GoogleAuthMode, referenceDiskLocalizationManifestFiles: Option[List[ValidFullGcsPath]]): PipelinesApiReferenceFilesMapping = {
    new PipelinesApiReferenceFilesMapping(auth, referenceDiskLocalizationManifestFiles)
  }

  def empty: PipelinesApiReferenceFilesMapping = new PipelinesApiReferenceFilesMapping(ApplicationDefaultMode("default"), None)
}
