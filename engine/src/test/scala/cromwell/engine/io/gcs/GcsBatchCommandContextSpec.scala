package cromwell.engine.io.gcs

import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import com.google.api.services.storage.StorageRequest
import common.assertion.CromwellTimeoutSpec
import cromwell.engine.io.gcs.GcsBatchCommandContextSpec.ExceptionSpewingGcsBatchIoCommand
import cromwell.filesystems.gcs.batch.GcsBatchIoCommand
import org.scalatest.BeforeAndAfter
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.util.{Failure, Success}

class GcsBatchCommandContextSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with Eventually with BeforeAndAfter  {
  behavior of "GcsBatchCommandContext"

  var commandContext: GcsBatchCommandContext[Unit, Unit] = _

  before {
    commandContext = GcsBatchCommandContext[Unit, Unit](
      request = ExceptionSpewingGcsBatchIoCommand(),
      replyTo = null
    )
  }

  it should "handle exceptions in success handlers" in {

    commandContext.promise.isCompleted should be(false)

    // Simulate a success response from an underlying IO operation:
    commandContext.callback.onSuccess((), new HttpHeaders())

    eventually {
      commandContext.promise.isCompleted should be(true)
    }

    commandContext.promise.future.value.get match {
      case Success(oops) => fail(s"Should not have produced a success: $oops")
      case Failure(error) => error.getMessage should be("Error processing IO response in onSuccessCallback: Mapping that value is impossible")
    }
  }

  it should "handle exceptions in failure handlers" in {

    commandContext.promise.isCompleted should be(false)

    // Simulate a success response from an underlying IO operation:
    commandContext.callback.onFailure(new GoogleJsonError { }, new HttpHeaders())

    eventually {
      commandContext.promise.isCompleted should be(true)
    }

    commandContext.promise.future.value.get match {
      case Success(oops) => fail(s"Should not have produced a success: $oops")
      case Failure(error) => error.getMessage should be("Error processing IO response in onFailureCallback: exception in failure handler...")
    }
  }
}

object GcsBatchCommandContextSpec {
  case class ExceptionSpewingGcsBatchIoCommand() extends GcsBatchIoCommand[Unit, Unit] {
    override def operation: StorageRequest[Unit] = ???
    override protected def mapGoogleResponse(response: Unit): Unit = throw new Exception("Mapping that value is impossible")
    override def withUserProject: GcsBatchIoCommand[Unit, Unit] = ???
    override def name: String = ???

    override def onFailure(googleJsonError: GoogleJsonError, httpHeaders: HttpHeaders): Option[Either[Unit, GcsBatchIoCommand[Unit, Unit]]] = {
      throw new Exception("exception in failure handler...")
    }

    override def commandDescription: String = s"ExceptionSpewingGcsBatchIoCommand (mock class for tests)"
  }
}
