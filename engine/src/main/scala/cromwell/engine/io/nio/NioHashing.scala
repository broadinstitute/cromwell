package cromwell.engine.io.nio

import cats.effect.IO
import com.google.cloud.storage.Blob
import common.util.StringUtil.EnhancedString
import cromwell.core.callcaching.HashType.HashType
import cromwell.core.callcaching.{FileHash, FileHashStrategy, HashType}
import cromwell.core.path.Path
import cromwell.filesystems.blob.BlobPath
import cromwell.filesystems.drs.DrsPath
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.http.HttpPath
import cromwell.filesystems.s3.S3Path
import cromwell.util.TryWithResource.tryWithResource

import scala.util.Try

object NioHashing {

  def hash(file: Path, hashStrategy: FileHashStrategy): IO[String] =
    // If there is no hash accessible from the file storage system,
    // we'll read the file and generate the hash ourselves if we can.
    getStoredHash(file, hashStrategy)
      .flatMap {
        case Some(storedHash) => IO.pure(storedHash)
        case None =>
          if (canHashLocally(file))
            generateMd5FileHashForPath(file)
          else
            IO.raiseError(
              new Exception(
                s"File of type ${file.getClass.getSimpleName} requires an associated hash, not present for ${file.pathAsString.maskSensitiveUri}"
              )
            )
      }
      .map(_.hash)

  def getStoredHash(file: Path, hashStrategy: FileHashStrategy): IO[Option[FileHash]] =
    file match {
      case gcsPath: GcsPath => getFileHashForGcsPath(gcsPath, hashStrategy)
      case blobPath: BlobPath => getFileHashForBlobPath(blobPath, hashStrategy)
      case drsPath: DrsPath => getFileHashForDrsPath(drsPath, hashStrategy)
      case s3Path: S3Path => getFileHashForS3Path(s3Path, hashStrategy)
      case _ => IO.pure(None)
    }

  /**
    * In some scenarios like SFS it is appropriate for Cromwell to hash files using its own CPU power.
    *
    * In cloud scenarios, we don't want this because the files are huge, downloading them is slow & expensive,
    * and the extreme CPU usage destabilizes the instance (WX-1566). For more context, see also comments
    * on `cromwell.filesystems.blob.BlobPath#largeBlobFileMetadataKey()`.
    *
    * Cromwell is fundamentally supposed to be a job scheduler, and heavy computation should take place elsewhere.
    *
    * @param file The path to consider for local hashing
    */
  private def canHashLocally(file: Path) =
    file match {
      case _: HttpPath => false
      case _: GcsPath => false
      case _: BlobPath => false
      case _: DrsPath => false
      case _: S3Path => false
      case _ => true
    }

  private def generateMd5FileHashForPath(path: Path): IO[FileHash] = delayedIoFromTry {
    tryWithResource(() => path.newInputStream) { inputStream =>
      FileHash(HashType.Md5, org.apache.commons.codec.digest.DigestUtils.md5Hex(inputStream))
    }
  }

  private def getFileHashForGcsPath(gcsPath: GcsPath, hashStrategy: FileHashStrategy): IO[Option[FileHash]] =
    delayedIoFromTry {
      val cloudFile = gcsPath.objectBlobId.map(id => gcsPath.cloudStorage.get(id))
      cloudFile.map(file =>
        hashStrategy.getFileHash(
          file,
          (f: Blob, hashType: HashType) =>
            hashType match {
              case HashType.Crc32c => Option(f.getCrc32c)
              case HashType.Md5 => Option(f.getMd5)
              case HashType.Identity => GcsPath.getBlobFingerprint(f)
              case _ => None
            }
        )
      )
    }

  private def getFileHashForBlobPath(blobPath: BlobPath, hashStrategy: FileHashStrategy): IO[Option[FileHash]] =
    IO {
      hashStrategy.getFileHash(
        blobPath,
        (f: BlobPath, hashType: HashType) =>
          hashType match {
            case HashType.Md5 => blobPath.md5HexString.toOption.flatten
            case _ => None
          }
      )
    }

  private def getFileHashForDrsPath(drsPath: DrsPath, hashStrategy: FileHashStrategy): IO[Option[FileHash]] =
    IO {
      val drsHashes = drsPath.getFileHashes
      hashStrategy.getFileHash(
        drsPath,
        (f: DrsPath, hashType: HashType) => drsHashes.get(hashType.toString.toLowerCase)
      )
    }

  private def getFileHashForS3Path(s3Path: S3Path, hashStrategy: FileHashStrategy): IO[Option[FileHash]] =
    IO {
      hashStrategy.getFileHash(
        s3Path,
        (f: S3Path, hashType: HashType) =>
          hashType match {
            case HashType.Etag => Option(s3Path.eTag)
            case _ => None
          }
      )
    }

  /**
    * Lazy evaluation of a Try in a delayed IO. This avoids accidentally eagerly evaluating the Try.
    *
    * IMPORTANT: Use this instead of IO.fromTry to make sure the Try will be reevaluated if the
    * IoCommand is retried.
    */
  private def delayedIoFromTry[A](t: => Try[A]): IO[A] = IO[A](t.get)
}
