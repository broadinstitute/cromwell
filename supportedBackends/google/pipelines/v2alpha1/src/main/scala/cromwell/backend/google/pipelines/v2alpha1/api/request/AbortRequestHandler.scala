package cromwell.backend.google.pipelines.v2alpha1.api.request

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.backend.google.pipelines.common.api.clients.PipelinesApiAbortClient.{PAPIAbortRequestSuccessful, PAPIOperationHasAlreadyFinished}
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

trait AbortRequestHandler extends LazyLogging { this: RequestHandler =>
  protected def handleGoogleError(abortQuery: PAPIAbortRequest, pollingManager: ActorRef, e: GoogleJsonError, responseHeaders: HttpHeaders): Try[Unit] = {
    // No need to fail the request if the job already terminated, we're all good
    // As seen in `cromwell.backend.google.pipelines.common.api.clients.PipelinesApiAbortClient` the difference between "finished" and "cancelled" only affects logging
    // Enhance to be more specific if/when Google implements https://partnerissuetracker.corp.google.com/issues/171993833
    if (Option(e.getCode).contains(400)) {
      logger.info(s"PAPI declined to abort job ${abortQuery.jobId.jobId} in workflow ${abortQuery.workflowId}, most likely because it is no longer running. Marking as finished. Message: ${e.getMessage}")
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
