package cromwell.backend.google.pipelines.v2alpha1.api.request

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.backend.google.pipelines.common.api.clients.PipelinesApiAbortClient.{PAPIAbortRequestSuccessful, PAPIOperationAlreadyCancelled, PAPIOperationHasAlreadyFinished}
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

trait AbortRequestHandler { this: RequestHandler =>
  protected def handleGoogleError(abortQuery: PAPIAbortRequest, pollingManager: ActorRef, e: GoogleJsonError, responseHeaders: HttpHeaders): Try[Unit] = {
    // No need to fail the request if the job was already cancelled, we're all good
    if (Option(e.getCode).contains(400) && Option(e.getMessage).contains("Operation has already been canceled")) {
      abortQuery.requester ! PAPIOperationAlreadyCancelled(abortQuery.jobId.jobId)
      Success(())
    } else if (Option(e.getCode).contains(400) && Option(e.getMessage).contains("Operation has already finished")) {
      abortQuery.requester ! PAPIOperationHasAlreadyFinished(abortQuery.jobId.jobId)
      Success(())
    } else {
      pollingManager ! PipelinesApiAbortQueryFailed(abortQuery, new SystemPAPIApiException(GoogleJsonException(e, responseHeaders)))
      Failure(new Exception(mkErrorString(e)))
    }
  }

  // The Genomics batch endpoint doesn't seem to be able to handle abort requests on V2 operations at the moment
  // For now, don't batch the request and execute it on its own 
  def handleRequest(abortQuery: PAPIAbortRequest, batch: BatchRequest, pollingManager: ActorRef)(implicit ec: ExecutionContext): Future[Try[Unit]] = {
    Future(abortQuery.httpRequest.execute()) map {
      case response if response.isSuccessStatusCode =>
        abortQuery.requester ! PAPIAbortRequestSuccessful(abortQuery.jobId.jobId)
        Success(())
      case response => for {
        asGoogleError <- Try(GoogleJsonError.parse(GoogleAuthMode.jsonFactory, response))
        handled <- handleGoogleError(abortQuery, pollingManager, asGoogleError, response.getHeaders)
      } yield handled
    } recover {
      case e =>
        pollingManager ! PipelinesApiAbortQueryFailed(abortQuery, new SystemPAPIApiException(e))
        Failure(e)
    }
  }
}
