package cromwell.core.callcaching

import cromwell.core.callcaching.HashType.HashType

import java.nio.{ByteBuffer, ByteOrder}
import java.security.MessageDigest
import java.util.Base64
import java.util.zip.CRC32C

// File hashing strategies used by IoHashCommand, primarily when obtaining file hashes
// for call caching purposes.
case class FileHashStrategy(priorityHashList: List[HashType]) {
  override def toString = s"FileHashStrategy(${priorityHashList.map(_.toString).mkString(", ")})"

  def getFileHash[A](fileToHash: A, hashFunc: (A, HashType) => Option[String]): Option[FileHash] =
    priorityHashList.flatMap(ht => hashFunc(fileToHash, ht).map(FileHash(ht, _))).headOption

}

object FileHashStrategy {
  val Crc32c: FileHashStrategy = FileHashStrategy(List(HashType.Crc32c))
  val Md5: FileHashStrategy = FileHashStrategy(List(HashType.Md5))
  val ETag: FileHashStrategy = FileHashStrategy(List(HashType.Etag))
  val Drs: FileHashStrategy = FileHashStrategy(List(HashType.Crc32c, HashType.Md5, HashType.Sha256, HashType.Etag))

  // TODO alert about bad hashes from config
  // TODO validate fs type here?
  def of(hashes: List[String]): FileHashStrategy = FileHashStrategy(hashes.flatMap(HashType(_)))
}

object HashType extends Enumeration {
  type HashType = Value

  // crc32c as a hex string
  val Crc32c: HashType.Value = Value
  val Md5: HashType.Value = Value
  val Identity: HashType.Value = Value
  val Etag: HashType.Value = Value
  val Sha256: HashType.Value = Value

  def apply(s: String): Option[HashType] = values.find(_.toString.toLowerCase == s.toLowerCase)

  implicit class HashTypeValue(hashType: Value) {
    def calculateHash(s: String, hexCrc32c: Boolean = false): String = hashType match {
      case Crc32c if hexCrc32c =>
        val crc32c = new CRC32C()
        crc32c.update(s.getBytes)
        crc32c.getValue.toHexString
      case Crc32c =>
        val crc32c = new CRC32C()
        crc32c.update(s.getBytes)
        val byteBuffer = ByteBuffer.allocate(4).order(ByteOrder.BIG_ENDIAN)
        byteBuffer.putInt(crc32c.getValue.toInt)
        Base64.getEncoder.encodeToString(byteBuffer.array)
      case Md5 =>
        org.apache.commons.codec.digest.DigestUtils.md5Hex(s)
      case Identity => s
      case Etag =>
        val chunkSize = 8 * 1024 * 1024
        val numChunks = (s.length.toDouble / chunkSize).ceil.toInt
        val parts = s.getBytes.grouped(chunkSize).map(org.apache.commons.codec.digest.DigestUtils.md5Hex)
        numChunks match {
          case 1 => parts.next()
          case _ =>
            s"${org.apache.commons.codec.digest.DigestUtils.md5Hex(parts.mkString)}-${numChunks}"
        }
      case Sha256 =>
        MessageDigest.getInstance("SHA-256").digest(s.getBytes).map("%02x" format _).mkString
    }
  }
}

case class FileHash(hashType: HashType, hash: String) {

  // GCS uses base64-encoded crc32c hashes (8 chars), other systems use hex (16)
  def isHexCrc32c = hashType == HashType.Crc32c && hash.length == 16

  /*
   * Compute the hash of the input string using a method equivalent to the one that produced this hash.
   */
  def computeHashOf(value: String): String =
    hashType.calculateHash(value, isHexCrc32c)
}
