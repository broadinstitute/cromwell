package cromwell.backend.google.pipelines.batch

import akka.actor.{Actor, ActorLogging, Props}

final class GcpBatchBackendSingletonActor extends Actor with ActorLogging{
  override def receive: Receive = ???
}

object GcpBatchBackendSingletonActor {
  def props(): Props = Props(new GcpBatchBackendSingletonActor)
}
