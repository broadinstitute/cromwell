package cromwell.engine.io.nio

import cats.effect.IO
import cloud.nio.spi.{FileHash, HashType}
import common.util.StringUtil.EnhancedString
import cromwell.core.path.Path
import cromwell.filesystems.blob.BlobPath
import cromwell.filesystems.drs.DrsPath
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.http.HttpPath
import cromwell.filesystems.s3.S3Path
import cromwell.util.TryWithResource.tryWithResource

import scala.util.Try

object NioHashing {

  def hash(file: Path): IO[String] =
    // If there is no hash accessible from the file storage system,
    // we'll read the file and generate the hash ourselves if we can.
    getStoredHash(file)
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

  def getStoredHash(file: Path): IO[Option[FileHash]] =
    file match {
      case gcsPath: GcsPath => getFileHashForGcsPath(gcsPath).map(Option(_))
      case blobPath: BlobPath => getFileHashForBlobPath(blobPath)
      case drsPath: DrsPath =>
        IO {
          // We assume all DRS files have a stored hash; this will throw
          // if the file does not.
          drsPath.getFileHash
        }.map(Option(_))
      case s3Path: S3Path =>
        IO {
          Option(FileHash(HashType.S3Etag, s3Path.eTag))
        }
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
      case _: BlobPath => false
      case _ => true
    }

  private def generateMd5FileHashForPath(path: Path): IO[FileHash] = delayedIoFromTry {
    tryWithResource(() => path.newInputStream) { inputStream =>
      FileHash(HashType.Md5, org.apache.commons.codec.digest.DigestUtils.md5Hex(inputStream))
    }
  }

  private def getFileHashForGcsPath(gcsPath: GcsPath): IO[FileHash] = delayedIoFromTry {
    gcsPath.objectBlobId.map(id => FileHash(HashType.GcsCrc32c, gcsPath.cloudStorage.get(id).getCrc32c))
  }

  private def getFileHashForBlobPath(blobPath: BlobPath): IO[Option[FileHash]] = delayedIoFromTry {
    blobPath.md5HexString.map(md5 => md5.map(FileHash(HashType.Md5, _)))
  }

  /**
    * Lazy evaluation of a Try in a delayed IO. This avoids accidentally eagerly evaluating the Try.
    *
    * IMPORTANT: Use this instead of IO.fromTry to make sure the Try will be reevaluated if the
    * IoCommand is retried.
    */
  private def delayedIoFromTry[A](t: => Try[A]): IO[A] = IO[A](t.get)
}
