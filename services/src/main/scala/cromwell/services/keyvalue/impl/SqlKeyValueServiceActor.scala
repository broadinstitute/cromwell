package cromwell.services.keyvalue.impl

import akka.actor.{ActorRef, Props}
import com.typesafe.config.Config
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.keyvalue.KeyValueServiceActor

import scala.concurrent.duration._
object SqlKeyValueServiceActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = Props(SqlKeyValueServiceActor(serviceConfig, globalConfig, serviceRegistryActor)).withDispatcher(ServiceDispatcher)
  val WriteQueueThreshold = 10 * 1000
  val ReadQueueThreshold = 10 * 1000
}

final case class SqlKeyValueServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
  extends KeyValueServiceActor {

  override protected def kvReadActorProps = {
    SqlKeyValueReadActor.props(
      SqlKeyValueServiceActor.WriteQueueThreshold,
      serviceRegistryActor
    )
  }

  override protected def kvWriteActorProps = {
    SqlKeyValueWriteActor.props(
      SqlKeyValueServiceActor.WriteQueueThreshold,
      serviceRegistryActor,
      5.seconds,
      1000
    )
  }
}
