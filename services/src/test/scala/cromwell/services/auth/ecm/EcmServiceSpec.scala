package cromwell.services.auth.ecm

import akka.http.scaladsl.model.StatusCodes
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

class EcmServiceSpec extends AnyFlatSpec with Matchers with TableDrivenPropertyChecks {

  private val ecmService = new EcmService("https://mock-ecm-url.org")

  private val ecm400ErrorMsg = "No enum constant bio.terra.externalcreds.generated.model.Provider.MyOwnProvider"
  private val ecm404ErrorMsg =
    "No linked account found for user ID: 123 and provider: github. Please go to the Terra Profile page External Identities tab to link your account for this provider"

  private val testCases = Table(
    ("test name", "response status code", "response body string", "expected error message"),
    ("return custom 401 error when status code is 401",
     StatusCodes.Unauthorized,
     "<h2>could be anything</h2>",
     "Invalid or missing authentication credentials."
    ),
    ("extract message from valid ErrorReport JSON if status code is 400",
     StatusCodes.BadRequest,
     s"""{ "message" : "$ecm400ErrorMsg", "statusCode" : 400}""",
     ecm400ErrorMsg
    ),
    ("extract message from valid ErrorReport JSON if status code is 404",
     StatusCodes.NotFound,
     s"""{ "message" : "$ecm404ErrorMsg", "statusCode" : 404}""",
     ecm404ErrorMsg
    ),
    ("extract message from valid ErrorReport JSON if status code is 500",
     StatusCodes.InternalServerError,
     """{ "message" : "Internal error", "statusCode" : 500}""",
     "Internal error"
    ),
    ("return response error body if it fails to parse JSON",
     StatusCodes.InternalServerError,
     "Response error - not a JSON",
     "Response error - not a JSON"
    ),
    ("return response error body if JSON doesn't contain 'message' key",
     StatusCodes.BadRequest,
     """{"non-message-key" : "error message"}""",
     """{"non-message-key" : "error message"}"""
    )
  )

  behavior of "extractErrorMessage in EcmService"

  forAll(testCases) { (testName, statusCode, responseBodyAsStr, expectedErrorMsg) =>
    it should testName in {
      assert(ecmService.extractErrorMessage(statusCode, responseBodyAsStr) == expectedErrorMsg)
    }
  }
}
