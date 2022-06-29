package cloud.nio.spi

import java.nio.{ByteBuffer, ByteOrder}
import java.security.MessageDigest
import java.util.Base64
import java.util.zip.CRC32C

object HashType extends Enumeration {
  type HashType = Value

  // crc32c as a hex string
  val Crc32c: HashType.Value = Value
  // GCS crc32c, which is base64-encoded instead of a hex string
  val GcsCrc32c: HashType.Value = Value
  val Md5: HashType.Value = Value
  // AWS S3 etag
  val S3Etag: HashType.Value = Value
  val Sha256: HashType.Value = Value

  implicit class HashTypeValue(hashType: Value) {
    def calculateHash(s: String): String = hashType match {
      case Crc32c =>
        val crc32c = new CRC32C()
        crc32c.update(s.getBytes)
        crc32c.getValue.toHexString
      case GcsCrc32c =>
        val crc32c = new CRC32C()
        crc32c.update(s.getBytes)
        val byteBuffer = ByteBuffer.allocate(4).order(ByteOrder.BIG_ENDIAN)
        byteBuffer.putInt(crc32c.getValue.toInt)
        Base64.getEncoder.encodeToString(byteBuffer.array)
      case Md5 =>
        org.apache.commons.codec.digest.DigestUtils.md5Hex(s)
      case S3Etag =>
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
