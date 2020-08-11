package cromwell.backend.google.pipelines.common

import java.util

import _root_.io.circe.generic.auto._
import _root_.io.circe.parser._
import cats.effect.IO
import cats.implicits._
import com.google.api.services.storage.StorageScopes
import com.google.cloud.storage.{BlobId, Storage, StorageOptions}
import com.google.common.io.BaseEncoding
import com.google.common.primitives.Longs
import cromwell.backend.google.pipelines.common.io.PipelinesApiReferenceFilesDisk
import cromwell.filesystems.gcs.GcsPathBuilder
import cromwell.filesystems.gcs.GcsPathBuilder.{InvalidFullGcsPath, ValidFullGcsPath}
import com.google.cloud.storage.Storage.{BlobField, BlobGetOption}
import cromwell.backend.google.pipelines.common.errors.InvalidGcsPathsInManifestFileException
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import org.slf4j.{Logger, LoggerFactory}

/**
 * This file contains logic related to public Broad reference files stored in public GCS bucket.
 * During instance creation it will read and parse Broad reference disk manifest files from GCS and populate internal
 * map containing relations between reference file paths and reference disks. This may take significant time, so the
 * best way to do it is during Cromwell startup.
 *
 */
case class PipelinesApiReferenceFilesMapping(validReferenceFilesMap: Map[String, PipelinesApiReferenceFilesDisk])

object PipelinesApiReferenceFilesMapping extends PipelinesApiReferenceFilesMappingOperations

// functionality extracted into the trait for testing purposes
protected trait PipelinesApiReferenceFilesMappingOperations {
  private val logger: Logger = LoggerFactory.getLogger(getClass)

  case class ReferenceFile(path: String, crc32c: Long)
  case class ManifestFile(imageIdentifier: String, diskSizeGb: Int, files: List[ReferenceFile])

  def generateReferenceFilesMapping(auth: GoogleAuthMode,
                                    referenceDiskLocalizationManifestFiles: Option[List[ValidFullGcsPath]]): PipelinesApiReferenceFilesMapping = {

    val validReferenceFilesMapIO = referenceDiskLocalizationManifestFiles match {
      case None =>
        IO.pure(Map.empty[String, PipelinesApiReferenceFilesDisk])
      case Some(manifestFilesPaths) =>
        val gcsClient = StorageOptions
          .newBuilder()
          .setCredentials(auth.credentials(Set(StorageScopes.DEVSTORAGE_READ_ONLY)))
          .build
          .getService
        val manifestFilesIo = manifestFilesPaths traverse { manifestFilePath => readReferenceDiskManifestFileFromGCS(gcsClient, manifestFilePath) }
        manifestFilesIo flatMap { manifestFiles =>
          manifestFiles
            .traverse(manifestFile => getMapOfValidReferenceFilePathsToDisks(gcsClient, manifestFile))
            .map(_.flatten.toMap)
        }
    }
    PipelinesApiReferenceFilesMapping(validReferenceFilesMapIO.unsafeRunSync())
  }

  def getReferenceDisksToMount(pipelinesApiReferenceFilesMapping: PipelinesApiReferenceFilesMapping,
                               inputFilePaths: Set[String]): List[PipelinesApiReferenceFilesDisk] = {
    pipelinesApiReferenceFilesMapping.validReferenceFilesMap.filterKeys(key => inputFilePaths.contains(s"gs://$key")).values.toList.distinct
  }

  protected def readReferenceDiskManifestFileFromGCS(gcsClient: Storage, gcsPath: ValidFullGcsPath): IO[ManifestFile] = {
    val manifestFileBlobIo = IO { gcsClient.get(BlobId.of(gcsPath.bucket, gcsPath.path.substring(1))) }
    manifestFileBlobIo flatMap { manifestFileBlob =>
      val jsonStringIo = IO { manifestFileBlob.getContent().map(_.toChar).mkString }
      jsonStringIo.flatMap(jsonStr => IO.fromEither(decode[ManifestFile](jsonStr)))
    }
  }

  private def getReferenceFileToValidatedGcsPathMap(referenceFiles: Set[ReferenceFile]): IO[Map[ReferenceFile, ValidFullGcsPath]] = {
    val filesAndValidatedPaths = referenceFiles.map {
      referenceFile => (referenceFile, GcsPathBuilder.validateGcsPath(s"gs://${referenceFile.path}"))
    }.toMap

    val filesWithValidPaths = filesAndValidatedPaths.collect {
      case (referenceFile, validPath: ValidFullGcsPath) => (referenceFile, validPath)
    }
    val filesWithInvalidPaths = filesAndValidatedPaths.collect {
      case (referenceFile, invalidPath: InvalidFullGcsPath) => (referenceFile, invalidPath)
    }

    if (filesWithInvalidPaths.nonEmpty) {
      IO.raiseError(new InvalidGcsPathsInManifestFileException(filesWithInvalidPaths.keySet.map(_.path).toList))
    } else {
      IO.pure(filesWithValidPaths)
    }
  }

  protected def bulkValidateCrc32cs(gcsClient: Storage,
                          filesWithValidPathsIo: IO[Map[ReferenceFile, ValidFullGcsPath]]): IO[Map[ReferenceFile, Boolean]] = {
    val gcsBatch = gcsClient.batch()
    val filesAndBlobResultsIo = filesWithValidPathsIo.map { filesWithValidPaths =>
      val filesAndBlobResults = filesWithValidPaths map {
        case (referenceFile, ValidFullGcsPath(bucket, path)) =>
          val blobGetResult = gcsBatch.get(BlobId.of(bucket, path.substring(1)), BlobGetOption.fields(BlobField.CRC32C))
          (referenceFile, blobGetResult)
      }
      gcsBatch.submit()
      filesAndBlobResults
    }

    filesAndBlobResultsIo.map { filesAndBlobResults =>
      filesAndBlobResults map {
        case (referenceFile, blobGetResult) =>
          val crc32cFromManifest = BaseEncoding.base64.encode(
            // drop 4 leading bytes from Long crc32c value
            // https://stackoverflow.com/a/25111119/1794750
            util.Arrays.copyOfRange(Longs.toByteArray(referenceFile.crc32c), 4, 8)
          )

          (referenceFile, crc32cFromManifest === blobGetResult.get().getCrc32c)
      }
    }
  }

  private def getMapOfValidReferenceFilePathsToDisks(gcsClient: Storage, manifestFile: ManifestFile): IO[Map[String, PipelinesApiReferenceFilesDisk]] = {
    val refDisk = PipelinesApiReferenceFilesDisk(manifestFile.imageIdentifier, manifestFile.diskSizeGb)
    val allReferenceFilesFromManifestMap = manifestFile.files.map(refFile => (refFile, refDisk)).toMap
    val referenceFilesWithValidPathsIo = getReferenceFileToValidatedGcsPathMap(allReferenceFilesFromManifestMap.keySet)
    val filesWithValidatedCrc32csIo = bulkValidateCrc32cs(gcsClient, referenceFilesWithValidPathsIo)

    val validReferenceFilesFromManifestMapIo = filesWithValidatedCrc32csIo map {
      filesWithValidatedCrc32c => {
        allReferenceFilesFromManifestMap.filterKeys(key => filesWithValidatedCrc32c.getOrElse(key, false))
      }
    }

    validReferenceFilesFromManifestMapIo map { validReferenceFilesFromManifestMap =>
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
