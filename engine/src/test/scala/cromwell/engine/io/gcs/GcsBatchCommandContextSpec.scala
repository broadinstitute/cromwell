package cromwell.engine.io.gcs

import akka.actor.ActorRef
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import common.assertion.CromwellTimeoutSpec
import cromwell.filesystems.gcs.batch._
import org.scalatest.BeforeAndAfter
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.util.{Failure, Success}

class GcsBatchCommandContextSpec
  extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with Eventually with BeforeAndAfter {
  behavior of "GcsBatchCommandContext"

  it should "handle exceptions in success handlers" in {
    val exceptionSpewingCommandContext = GcsBatchCommandContext[Unit, Void](
      request = ExceptionSpewingGcsBatchIoCommand,
      replyTo = ActorRef.noSender
    )

    exceptionSpewingCommandContext.promise.isCompleted should be(false)

    // Simulate a success response from an underlying IO operation:
    exceptionSpewingCommandContext.callback.onSuccess(null, new HttpHeaders())

    eventually {
      exceptionSpewingCommandContext.promise.isCompleted should be(true)
    }

    exceptionSpewingCommandContext.promise.future.value.get match {
      case Success(oops) => fail(s"Should not have produced a success: $oops")
      case Failure(error) => error.getMessage should be("Error processing IO response in onSuccessCallback: Ill behaved code that throws in mapGoogleResponse")
    }
  }

  it should "handle exceptions in failure handlers" in {
    val exceptionSpewingCommandContext = GcsBatchCommandContext[Unit, Void](
      request = ExceptionSpewingGcsBatchIoCommand,
      replyTo = ActorRef.noSender
    )

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
    val errorReturningCommandContext = GcsBatchCommandContext[Unit, Unit](
      request = ErrorReturningGcsBatchIoCommand(),
      replyTo = ActorRef.noSender
    )

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
