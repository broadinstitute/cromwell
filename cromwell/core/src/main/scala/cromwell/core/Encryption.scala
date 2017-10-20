package cromwell.core

import java.security.SecureRandom
import javax.crypto.Cipher
import javax.crypto.spec.{IvParameterSpec, SecretKeySpec}

import org.apache.commons.codec.binary.Base64

import scala.util.{Failure, Success, Try}

/**
 * 256-bit AES encryption class for encryption and decryption.
 *
 * 'iv' is short for 'initialization vector'
 */
case object Aes256Cbc {

  val encryption = "AES"
  val blockSize = 128
  val keySize = 256
  val paddingMode = "PKCS5Padding"
  val cipherMode = "CBC"

  val ranGen = new SecureRandom()

  final def init(mode: Int, secretKey: Array[Byte], iv: Array[Byte]) = {
    val key = new SecretKeySpec(secretKey, encryption)
    val ivSpec = new IvParameterSpec(iv)
    val cipher = Cipher.getInstance(s"$encryption/$cipherMode/$paddingMode")
    cipher.init(mode, key, ivSpec)
    cipher
  }

  final def validateLength(arrayName: String, array: Array[Byte], expectedBitLength: Int): Try[Unit] = {
    if (array.length * 8 == expectedBitLength) {
      Success(())
    } else {
      Failure(new IllegalArgumentException(s"$arrayName size (${array.length * 8} bits) did not match the required length $expectedBitLength"))
    }
  }

  final def encrypt(plainText: Array[Byte], secretKey: SecretKey): Try[EncryptedBytes] = {
    validateLength("Secret key", secretKey.key, keySize) map { _ =>
      val iv = new Array[Byte](blockSize / 8)
      ranGen.nextBytes(iv)

      val cipher = init(Cipher.ENCRYPT_MODE, secretKey.key, iv)
      EncryptedBytes(cipher.doFinal(plainText), iv)
    }
  }

  final def decrypt(encryptedBytes: EncryptedBytes, secretKey: SecretKey): Try[Array[Byte]] = {
    for {
      _ <- validateLength("Secret key", secretKey.key, keySize)
      _ <- validateLength("Initialization vector", encryptedBytes.initializationVector, blockSize)
      bytes = init(Cipher.DECRYPT_MODE, secretKey.key, encryptedBytes.initializationVector).doFinal(encryptedBytes.cipherText)
    } yield bytes
  }
}

final case class EncryptedBytes(cipherText: Array[Byte], initializationVector: Array[Byte]) {
  def base64CipherText = Base64.encodeBase64String(cipherText)
  def base64Iv = Base64.encodeBase64String(initializationVector)
}

object EncryptedBytes {
  def apply(base64CipherTextString: String, base64IvString: String): EncryptedBytes =
    EncryptedBytes(Base64.decodeBase64(base64CipherTextString), Base64.decodeBase64(base64IvString))
}

final case class SecretKey(key: Array[Byte])
object SecretKey {
  def apply(base64KeyString: String): SecretKey = SecretKey(Base64.decodeBase64(base64KeyString))
}