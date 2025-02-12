package cromwell.core.callcaching

import com.typesafe.scalalogging.LazyLogging
import cromwell.core.callcaching.HashType.HashType

// File hashing strategies used by IoHashCommand, primarily when obtaining file hashes for call caching purposes.
// When obtaining the hash for a file we try each of the hash types in priorityHashList, in the order they appear
// in the list.
case class FileHashStrategy(priorityHashList: List[HashType]) {
  override def toString = s"FileHashStrategy(${priorityHashList.map(_.toString).mkString(", ")})"

  def isEmpty: Boolean = priorityHashList.isEmpty

  // Lazily evaluate hashes from `priorityList` until we find one that exists
  def getFileHash[A](fileToHash: A, hashFunc: (A, HashType) => Option[String]): Option[FileHash] =
    priorityHashList.to(LazyList).flatMap(ht => hashFunc(fileToHash, ht).map(FileHash(ht, _))).headOption
}

object FileHashStrategy extends LazyLogging {
  // Define some commonly-used strategies so they're easier to refer to
  val Crc32c: FileHashStrategy = FileHashStrategy(List(HashType.Crc32c))
  val Md5: FileHashStrategy = FileHashStrategy(List(HashType.Md5))
  val ETag: FileHashStrategy = FileHashStrategy(List(HashType.Etag))
  val Drs: FileHashStrategy = FileHashStrategy(List(HashType.Crc32c, HashType.Md5, HashType.Sha256, HashType.Etag))
  val Empty: FileHashStrategy = FileHashStrategy(List.empty)

  // No effort is made here to validate which has strategies are supported by which filesystems.
  // Developers who want to see the code that actually computes hashes should look in:
  // * cromwell.engine.io.nio.NioHashing
  // * cromwell.filesystems.gcs.batch.GcsBatchIoCommand
  // * cromwell.filesystems.s3.batch.S3BatchIoCommand
  def of(hashTypes: List[String]): FileHashStrategy = {
    val typedHashTypes = hashTypes.flatMap { h =>
      HashType(h) match {
        case Some(t) => Some(t)
        case None => logger.warn(s"Found invalid call caching hash type: $h"); None
      }
    }
    FileHashStrategy(typedHashTypes)
  }
}

object HashType extends Enumeration {
  type HashType = Value

  val Crc32c: HashType.Value = Value
  val Md5: HashType.Value = Value
  val Identity: HashType.Value = Value
  val Etag: HashType.Value = Value
  val Sha256: HashType.Value = Value

  def apply(s: String): Option[HashType] = values.find(_.toString.toLowerCase == s.toLowerCase)
}

case class FileHash(hashType: HashType, hash: String)
