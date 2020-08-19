package cromwell.cloudsupport.gcp.auth

import java.io.StringWriter
import java.security.KeyPairGenerator
import java.util.Base64

import com.google.api.client.http.{HttpHeaders, HttpResponseException}
import com.google.auth.oauth2.GoogleCredentials
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode.OptionLookup
import org.scalatest.Assertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

import scala.util.{Failure, Try}


class GoogleAuthModeSpec extends AnyFlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "GoogleAuthMode"

  private def mockHttpResponseException(statusCode: Int): HttpResponseException = {
    new HttpResponseException.Builder(statusCode, "mock message", new HttpHeaders).build()
  }

  private val testedExceptions = Table(
    ("description", "exception", "isFatal"),
    ("HTTP 400", mockHttpResponseException(400), true),
    ("HTTP 401", mockHttpResponseException(401), true),
    ("HTTP 402", mockHttpResponseException(402), false),
    ("HTTP 403", mockHttpResponseException(403), true),
    ("a random exception", new IllegalArgumentException(), false)
  )

  forAll(testedExceptions) { (description, exception, isFatal) =>
    it should s"return isFatal == $isFatal for $description" in {
      GoogleAuthMode.isFatal(exception) should be(isFatal)
    }
  }

  it should "create SimpleClientSecrets" in {
    val secrets = SimpleClientSecrets("id", "secret")
    secrets.clientId should be("id")
    secrets.clientSecret should be("secret")
  }
}

object GoogleAuthModeSpec {
  def assumeHasApplicationDefaultCredentials(): Unit = {
    tryApplicationDefaultCredentials match {
      case Failure(exception) => cancel(exception.getMessage)
      case _ =>
    }
    ()
  }

  private lazy val tryApplicationDefaultCredentials: Try[Unit] = Try {
    GoogleCredentials.getApplicationDefault
    ()
  }

  private def toJson(contents: (String, String)*): String = {
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

  lazy val userCredentialsContents: String = {
    toJson(
      "type" -> "authorized_user",
      "client_id" -> "the_id",
      "client_secret" -> "the_secret",
      "refresh_token" -> "the_token"
    )
  }

  lazy val serviceAccountPemContents: String = {
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

  lazy val serviceAccountJsonContents: String = {
    toJson(
      "type" -> "service_account",
      "client_id" -> "the_account_id",
      "client_email" -> "the_email",
      theStringThatShallNotBeNamed -> serviceAccountPemContents,
      s"${theStringThatShallNotBeNamed}_id" -> "the_key_id"
    )
  }

  lazy val refreshTokenOptions: OptionLookup = Map("refresh_token" -> "the_refresh_token")
  lazy val userServiceAccountOptions: OptionLookup = Map("user_service_account_json" -> serviceAccountJsonContents)
}
