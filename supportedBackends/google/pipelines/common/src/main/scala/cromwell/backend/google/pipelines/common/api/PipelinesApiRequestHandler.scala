package cromwell.backend.google.pipelines.common.api

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.PAPIApiRequest

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

trait PipelinesApiRequestHandler {
  def makeBatchRequest: BatchRequest
  def enqueue[T <: PAPIApiRequest](papiApiRequest: T, batchRequest: BatchRequest, pollingManager: ActorRef)(implicit ec: ExecutionContext): Future[Try[Unit]]
}
