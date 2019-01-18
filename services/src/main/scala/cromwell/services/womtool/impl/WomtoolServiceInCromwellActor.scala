package cromwell.services.womtool.impl

import akka.actor.{Actor, ActorRef, Props}
import akka.pattern.pipe
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.womtool.Describer
import cromwell.services.womtool.WomtoolServiceMessages._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.{ExecutionContext, Future}

class WomtoolServiceInCromwellActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with LazyLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive: Receive = {
    case DescribeRequest(filesCollection) =>

      // We are consciously wrapping a Future around the Await.result way down in the HTTP import resolver until we can update the whole call hierarchy to async
      // https://doc.akka.io/docs/akka/2.5.16/actors.html?language=scala#ask-send-and-receive-future
      Future {
        Describer.describeWorkflow(filesCollection)
      } pipeTo sender()
      ()
    case ShutdownCommand =>
      // This service doesn't do any work on shutdown but the service registry pattern requires it (see #2575)
      context.stop(self)
  }

}

object WomtoolServiceInCromwellActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) =
    Props(new WomtoolServiceInCromwellActor(serviceConfig, globalConfig, serviceRegistryActor)).withDispatcher(ServiceDispatcher)
}
