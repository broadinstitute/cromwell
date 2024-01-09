package cromwell.cloudsupport.gcp.auth

import com.google.api.client.http.{HttpHeaders, HttpResponseException}
import common.assertion.CromwellTimeoutSpec
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode.OptionLookup
import org.scalatest.Assertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

import scala.util.{Failure, Try}

class GoogleAuthModeSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "GoogleAuthMode"

  private def mockHttpResponseException(statusCode: Int): HttpResponseException =
    new HttpResponseException.Builder(statusCode, "mock message", new HttpHeaders).build()

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

object GoogleAuthModeSpec extends ServiceAccountTestSupport {
  def assumeHasApplicationDefaultCredentials(): Unit = {
    tryApplicationDefaultCredentials match {
      case Failure(exception) => cancel(exception.getMessage)
      case _ =>
    }
    ()
  }

  private lazy val tryApplicationDefaultCredentials: Try[Unit] = Try {
    new ApplicationDefaultMode("application-default").credentials()
    ()
  }

  lazy val userCredentialsContents: String =
    toJson(
      "type" -> "authorized_user",
      "client_id" -> "the_id",
      "client_secret" -> "the_secret",
      "refresh_token" -> "the_token"
    )

  lazy val refreshTokenOptions: OptionLookup = Map("refresh_token" -> "the_refresh_token")
  lazy val userServiceAccountOptions: OptionLookup = Map("user_service_account_json" -> serviceAccountJsonContents)
}
