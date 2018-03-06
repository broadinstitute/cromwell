package cromwell.services.loadcontroller.impl

import akka.actor.{Actor, ActorLogging, ActorRef, Timers}
import akka.routing.Listeners
import cats.data.NonEmptyList
import com.typesafe.config.Config
import cromwell.core.actor.BatchActor.QueueWeight
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.loadcontroller.LoadControllerService._
import cromwell.services.loadcontroller.LoadMetric
import cromwell.services.loadcontroller.impl.LoadControllerServiceActor._
import cromwell.services.metadata.MetadataService.ListenToMetadataWriteActor
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._

object LoadControllerServiceActor {
  val LoadControllerServiceName = "LoadController"
  val LoadLevelInstrumentation = NonEmptyList.one("loadLevel")
  case object LoadControlTimerKey
  case object LoadControlTimerAction
}

/**
  * Service Actor that monitors load and sends alert to the rest of the system when load is determined abnormal.
  */
class LoadControllerServiceActor(serviceConfig: Config, globalConfig: Config, override val serviceRegistryActor: ActorRef) extends Actor
  with ActorLogging with Listeners with Timers with CromwellInstrumentation {
  private val controlFrequency = serviceConfig.as[Option[FiniteDuration]]("control-frequency").getOrElse(5.seconds)
  private val metadataQueueThreshold = serviceConfig.as[Option[Int]]("metadata-queue-threshold").getOrElse(1 * 1000 * 1000)
  private val metadataQueueMetric = MetadataQueueMetric(metadataQueueThreshold)

  private var loadLevel: LoadLevel = NormalLoad
  private var loadMetrics: Map[LoadMetric, LoadLevel] = Map.empty

  override def receive = listenerManagement.orElse(controlReceive)

  override def preStart() = {
    serviceRegistryActor ! ListenToMetadataWriteActor
    timers.startPeriodicTimer(LoadControlTimerKey, LoadControlTimerAction, controlFrequency)
    super.preStart()
  }

  private def controlReceive: Receive = {
    case LoadControlTimerAction => checkLoad()
    case QueueWeight(weight) => updateMetric(metadataQueueMetric, weight)
    case ShutdownCommand => context stop self
  }

  def updateMetric(metric: LoadMetric, loadLevel: Int) = {
    loadMetrics = loadMetrics + (metric -> metric.loadLevel(loadLevel))
    sendGauge(NonEmptyList.of("load", metric.name), loadLevel.toLong, InstrumentationPrefixes.ServicesPrefix)
  }

  def checkLoad(): Unit = {
    // Simply take the max level of all load metrics for now
    val newLoadLevel = if (loadMetrics.nonEmpty) loadMetrics.values.max else NormalLoad
    // The load level escalates if the new load is higher than the previous load
    val escalates = loadLevelOrdering.gt(newLoadLevel, loadLevel)
    // Back to normal if we were not normal before but now are
    val backToNormal = loadLevel != NormalLoad && newLoadLevel == NormalLoad
    // If there's something to say, let it out !
    if (escalates || backToNormal) gossip(newLoadLevel)
    loadLevel = newLoadLevel
    sendGauge(LoadLevelInstrumentation, loadLevel.level.toLong)
  }
}
