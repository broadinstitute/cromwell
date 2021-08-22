package cromwell.backend.google.pipelines.common.api

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.http.HttpRequest
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.BatchRequestTimeoutConfiguration
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.PAPIApiRequest

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

trait PipelinesApiRequestHandler {
  def initializeHttpRequest(batchRequestTimeoutConfiguration: BatchRequestTimeoutConfiguration)
                           (httpRequest: HttpRequest): Unit = {
    batchRequestTimeoutConfiguration.readTimeoutMillis foreach {
      timeout => httpRequest.setReadTimeout(timeout.value)
    }
    batchRequestTimeoutConfiguration.connectTimeoutMillis foreach {
      timeout => httpRequest.setConnectTimeout(timeout.value)
    }
  }

  def makeBatchRequest: BatchRequest
  def enqueue[T <: PAPIApiRequest](papiApiRequest: T, batchRequest: BatchRequest, pollingManager: ActorRef)(implicit ec: ExecutionContext): Future[Try[Unit]]
}
