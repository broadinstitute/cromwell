package cromwell.engine.io.gcs

import cats.syntax.validated._
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import com.google.api.services.storage.StorageRequest
import cromwell.engine.io.gcs.GcsBatchCommandContextSpec.ExceptionSpewingGcsBatchIoCommand
import common.assertion.CromwellTimeoutSpec
import common.validation.ErrorOr.ErrorOr
import cromwell.engine.io.gcs.GcsBatchCommandContextSpec.{ErrorReturningGcsBatchIoCommand, ExceptionSpewingGcsBatchIoCommand}
import cromwell.filesystems.gcs.batch.GcsBatchIoCommand
import org.scalatest.BeforeAndAfter
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.util.{Failure, Success}

class GcsBatchCommandContextSpec extends AnyFlatSpec with Matchers with Eventually with BeforeAndAfter  {
  behavior of "GcsBatchCommandContext"

  var exceptionSpewingCommandContext: GcsBatchCommandContext[Unit, Unit] = _

  var errorReturningCommandContext: GcsBatchCommandContext[Unit, Unit] = _

  before {
    exceptionSpewingCommandContext = GcsBatchCommandContext[Unit, Unit](
      request = ExceptionSpewingGcsBatchIoCommand(),
      replyTo = null
    )

    errorReturningCommandContext = GcsBatchCommandContext[Unit, Unit](
      request = ErrorReturningGcsBatchIoCommand(),
      replyTo = null
    )
  }

  it should "handle exceptions in success handlers" in {

    exceptionSpewingCommandContext.promise.isCompleted should be(false)

    // Simulate a success response from an underlying IO operation:
    exceptionSpewingCommandContext.callback.onSuccess((), new HttpHeaders())

    eventually {
      exceptionSpewingCommandContext.promise.isCompleted should be(true)
    }

    exceptionSpewingCommandContext.promise.future.value.get match {
      case Success(oops) => fail(s"Should not have produced a success: $oops")
      case Failure(error) => error.getMessage should be("Error processing IO response in onSuccessCallback: Ill behaved code that throws in mapGoogleResponse")
    }
  }

  it should "handle exceptions in failure handlers" in {

    exceptionSpewingCommandContext.promise.isCompleted should be(false)

    // Simulate a failure response from an underlying IO operation:
    exceptionSpewingCommandContext.callback.onFailure(new GoogleJsonError { }, new HttpHeaders())

    eventually {
      exceptionSpewingCommandContext.promise.isCompleted should be(true)
    }

    exceptionSpewingCommandContext.promise.future.value.get match {
      case Success(oops) => fail(s"Should not have produced a success: $oops")
      case Failure(error) => error.getMessage should be("Error processing IO response in onFailureCallback: Ill behaved code that throws in onFailure")
    }
  }

  it should "handle errors in onSuccess" in {
    errorReturningCommandContext.promise.isCompleted should be(false)

    // Simulate a success response from an underlying IO operation:
    errorReturningCommandContext.callback.onSuccess((), new HttpHeaders())

    eventually {
      errorReturningCommandContext.promise.isCompleted should be(true)
    }

    errorReturningCommandContext.promise.future.value.get match {
      case Success(oops) => fail(s"Should not have produced a success: $oops")
      case Failure(error) => error.getMessage should be("Unexpected result in successful Google API call:\nWell behaved code that returns an error in mapGoogleResponse")
    }
  }
}

object GcsBatchCommandContextSpec {
  case class ExceptionSpewingGcsBatchIoCommand() extends GcsBatchIoCommand[Unit, Unit] {
    override def operation: StorageRequest[Unit] = ???
    override protected def mapGoogleResponse(response: Unit): ErrorOr[Unit] = throw new Exception("Ill behaved code that throws in mapGoogleResponse")
    override def withUserProject: GcsBatchIoCommand[Unit, Unit] = ???
    override def name: String = ???

    override def onFailure(googleJsonError: GoogleJsonError, httpHeaders: HttpHeaders): Option[Either[Unit, GcsBatchIoCommand[Unit, Unit]]] = {
      throw new Exception("Ill behaved code that throws in onFailure")
    }

    override def commandDescription: String = s"ExceptionSpewingGcsBatchIoCommand (mock class for tests)"
  }

  case class ErrorReturningGcsBatchIoCommand() extends GcsBatchIoCommand[Unit, Unit] {
    override def operation: StorageRequest[Unit] = ???
    override protected def mapGoogleResponse(response: Unit): ErrorOr[Unit] = "Well behaved code that returns an error in mapGoogleResponse".invalidNel
    override def withUserProject: GcsBatchIoCommand[Unit, Unit] = ???
    override def name: String = ???

    override def onSuccess(response: Unit, httpHeaders: HttpHeaders): ErrorOr[Either[Unit, GcsBatchIoCommand[Unit, Unit]]] = super.onSuccess(response, httpHeaders)
  }
}
