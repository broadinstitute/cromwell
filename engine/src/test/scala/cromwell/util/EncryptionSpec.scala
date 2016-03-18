package cromwell.util

import javax.crypto.Cipher

import cromwell.core.{SecretKey, Aes256Cbc}
import org.scalatest.{Assertions, FlatSpec, Matchers}

object EncryptionSpec {
  import Assertions._
  def assumeAes256Cbc() = assume(
    Cipher.getMaxAllowedKeyLength(Aes256Cbc.encryption) >= Aes256Cbc.keySize,
    """|
       |Did you install the Java Cryptography Extension (JCE) Unlimited Strength Jurisdiction Policy Files?
       |http://www.oracle.com/technetwork/java/javase/downloads/jce8-download-2133166.html
       |""".stripMargin)
}

class EncryptionSpec extends FlatSpec with Matchers {
  "Aes256Cbc" should "encrypt and decrypt" in {
    EncryptionSpec.assumeAes256Cbc()
    val charSet = "utf-8"
    val plainText = "Hello world".getBytes(charSet)
    val keyBytes = new Array[Byte](Aes256Cbc.keySize / 8)
    Aes256Cbc.ranGen.nextBytes(keyBytes)
    val secretKey = new SecretKey(keyBytes)
    val tryEncryptDecrypt = for {
      encryptedBytes <- Aes256Cbc.encrypt(plainText, secretKey)
      decryptedBytes <- Aes256Cbc.decrypt(encryptedBytes, secretKey)
    } yield new String(decryptedBytes, charSet)
    tryEncryptDecrypt.get should be("Hello world")
  }
}
