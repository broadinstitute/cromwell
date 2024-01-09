package cromwell.services.keyvalue.impl

import akka.actor.{ActorRef, Props}
import com.typesafe.config.Config
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.LoadConfig
import cromwell.services.keyvalue.KeyValueServiceActor
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
object SqlKeyValueServiceActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) =
    Props(SqlKeyValueServiceActor(serviceConfig, globalConfig, serviceRegistryActor)).withDispatcher(ServiceDispatcher)
}

final case class SqlKeyValueServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
    extends KeyValueServiceActor {

  private lazy val dbFlushRate = serviceConfig.as[Option[FiniteDuration]]("db-flush-rate").getOrElse(5.seconds)
  private lazy val dbBatchSize = serviceConfig.as[Option[Int]]("db-batch-size").getOrElse(200)

  override protected def kvReadActorProps =
    SqlKeyValueReadActor.props(
      LoadConfig.KeyValueReadThreshold,
      serviceRegistryActor
    )

  override protected def kvWriteActorProps =
    SqlKeyValueWriteActor.props(
      LoadConfig.KeyValueWriteThreshold,
      serviceRegistryActor,
      dbFlushRate,
      dbBatchSize
    )
}
