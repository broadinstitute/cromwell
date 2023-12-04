package cromwell.backend.google.pipelines.common.api

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.PAPIApiRequest

import scala.concurrent.Future
import scala.util.Try

trait PipelinesApiBatchHandler {
  def makeBatchRequest: BatchRequest
  def enqueue[T <: PAPIApiRequest](papiApiRequest: T,
                                   batchRequest: BatchRequest,
                                   pollingManager: ActorRef
  ): Future[Try[Unit]]
}
