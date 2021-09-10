package cromwell.cloudsupport.gcp.auth

import java.io.StringWriter
import java.security.KeyPairGenerator
import java.util.Base64

trait ServiceAccountTestSupport {
  def serviceAccountPemContents: String = {
    val keyPairGenerator = KeyPairGenerator.getInstance("RSA")
    keyPairGenerator.initialize(1024)
    val keyPair = keyPairGenerator.genKeyPair

    // extract the encoded private key, this is an unencrypted PKCS#8 private key
    val privateKey = keyPair.getPrivate
    val byteEncoded = privateKey.getEncoded
    val base64Encoded = Base64.getEncoder.encodeToString(byteEncoded)
    s"""|-----BEGIN PRIVATE KEY-----
        |$base64Encoded
        |-----END PRIVATE KEY-----
        |""".stripMargin
  }

  // Hide me from git secrets false positives
  private val theStringThatShallNotBeNamed = List("private", "key").mkString("_")

  def serviceAccountJsonContents: String = {
    toJson(
      "type" -> "service_account",
      "client_id" -> "the_account_id",
      "client_email" -> "the_email",
      theStringThatShallNotBeNamed -> serviceAccountPemContents,
      s"${theStringThatShallNotBeNamed}_id" -> "the_key_id"
    )
  }

  def toJson(contents: (String, String)*): String = {
    // Generator doesn't matter as long as it generates JSON. Using `jsonFactory` to get an extra line hit of coverage.
    val factory = GoogleAuthMode.jsonFactory
    val writer = new StringWriter()
    val generator = factory.createJsonGenerator(writer)
    generator.enablePrettyPrint()
    generator.writeStartObject()
    contents foreach {
      case (key, value) =>
        generator.writeFieldName(key)
        generator.writeString(value)
    }
    generator.writeEndObject()
    generator.close()
    writer.toString
  }
}
