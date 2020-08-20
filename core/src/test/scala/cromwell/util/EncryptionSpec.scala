package cromwell.util

import cromwell.core.{Aes256Cbc, SecretKey}
import javax.crypto.Cipher
import org.scalatest.Assertions
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

object EncryptionSpec {
  import Assertions._
  def assumeAes256Cbc() = assume(
    Cipher.getMaxAllowedKeyLength(Aes256Cbc.encryption) >= Aes256Cbc.keySize,
    """|
       |Did you install the Java Cryptography Extension (JCE) Unlimited Strength Jurisdiction Policy Files?
       |http://www.oracle.com/technetwork/java/javase/downloads/jce8-download-2133166.html
       |""".stripMargin)
}

class EncryptionSpec extends AnyFlatSpec with Matchers {
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
