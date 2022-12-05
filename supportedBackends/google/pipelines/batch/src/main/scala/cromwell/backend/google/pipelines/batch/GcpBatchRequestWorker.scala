package cromwell.backend.google.pipelines.batch

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.services.instrumentation.CromwellInstrumentationActor

import scala.concurrent.ExecutionContext

class GcpBatchRequestWorker(override val serviceRegistryActor: ActorRef) extends Actor with ActorLogging with CromwellInstrumentationActor {

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive = {
    case PipelinesApiWorkBatch(workBatch) =>

  }

}
