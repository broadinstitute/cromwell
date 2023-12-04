package cromiam.server.status

import com.softwaremill.sttp.asynchttpclient.future.AsyncHttpClientFutureBackend
import com.softwaremill.sttp.{sttp, Uri}
import org.broadinstitute.dsde.workbench.util.health.SubsystemStatus

import scala.concurrent.{ExecutionContext, Future}

/**
  * Represents a workbench system in use by CaaS which is checked via a call to its "status" endpoint
  */
trait StatusCheckedSubsystem {
  val statusUri: Uri

  implicit val sttpBackend = AsyncHttpClientFutureBackend()

  /**
    * Make a call to the status endpoint. If we receive a 200 OK fill in the SubsystemStatus w/ OK = true and no
    * error messages, otherwise OK = false and include the response body
    */
  def subsystemStatus()(implicit ec: ExecutionContext): Future[SubsystemStatus] =
    sttp.get(statusUri).send() map { x =>
      x.body match {
        case Right(_) => SubsystemStatus(true, None)
        case Left(errors) => SubsystemStatus(false, Option(List(errors)))
      }
    }
}
