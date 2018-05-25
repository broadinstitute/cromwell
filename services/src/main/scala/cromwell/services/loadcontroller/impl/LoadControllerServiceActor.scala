package cromwell.services.loadcontroller.impl

import akka.actor.{Actor, ActorLogging, ActorRef, Terminated, Timers}
import akka.routing.Listeners
import cats.data.NonEmptyList
import com.typesafe.config.Config
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.loadcontroller.LoadControllerService._
import cromwell.services.loadcontroller.impl.LoadControllerServiceActor._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._

object LoadControllerServiceActor {
  val LoadControllerServiceName = "LoadController"
  val LoadInstrumentationPrefix = Option("load")
  case object LoadControlTimerKey
  case object LoadControlTimerAction
  case class ActorAndMetric(actorRef: ActorRef, metricPath: NonEmptyList[String])
}

/**
  * Actor that monitors load and sends alert to the rest of the system when load is determined abnormal.
  * Send a LoadMetric message to the serviceRegistryActor to let this actor know about your current load.
  * Send a ListenToLoadController message to the serviceRegistryActor to listen to load updates.
  */
class LoadControllerServiceActor(serviceConfig: Config,
                                 globalConfig: Config,
                                 override val serviceRegistryActor: ActorRef
                                ) extends Actor
  with ActorLogging with Listeners with Timers with CromwellInstrumentation {
  private val controlFrequency = serviceConfig
    .as[Option[Duration]]("control-frequency")
    .getOrElse(5.seconds)

  private [impl] var loadLevel: LoadLevel = NormalLoad
  private [impl] var monitoredActors: Set[ActorRef] = Set.empty
  private [impl] var loadMetrics: Map[ActorAndMetric, LoadLevel] = Map.empty
  
  override def receive = listenerManagement.orElse(controlReceive)

  override def preStart() = {
    if (controlFrequency.isFinite()) 
      timers.startPeriodicTimer(LoadControlTimerKey, LoadControlTimerAction, controlFrequency.asInstanceOf[FiniteDuration])
    else 
      log.info("Load control disabled")
    super.preStart()
  }

  private val controlReceive: Receive = {
    case metric: LoadMetric => updateMetric(metric)
    case LoadControlTimerAction => checkLoad()
    case Terminated(terminee) => handleTerminated(terminee)
    case ShutdownCommand => context stop self
  }

  private def updateMetric(metric: LoadMetric): Unit = {
    val snd = sender()
    loadMetrics = loadMetrics + (ActorAndMetric(snd, metric.name) -> metric.loadLevel)
    if (!monitoredActors.contains(snd)) {
      context.watch(snd)
      monitoredActors = monitoredActors + snd
    }
    sendGauge(metric.name, metric.loadLevel.level.toLong, LoadInstrumentationPrefix)
  }

  private def checkLoad(): Unit = {
    // Simply take the max level of all load metrics for now
    val newLoadLevel = if (loadMetrics.nonEmpty) loadMetrics.values.max else NormalLoad
    // The load level escalates if the new load is higher than the previous load
    val escalates = loadLevelOrdering.gt(newLoadLevel, loadLevel)
    // Back to normal if we were not normal before but now are
    val backToNormal = loadLevel != NormalLoad && newLoadLevel == NormalLoad
    // If there's something to say, let it out !
    if (escalates || backToNormal) {
      if (newLoadLevel == HighLoad) log.info(s"The following components have reported being overloaded: $highLoadMetricsForLogging")
      gossip(newLoadLevel)
    }
    loadLevel = newLoadLevel
    sendGauge(NonEmptyList.one("global"), loadLevel.level.toLong, LoadInstrumentationPrefix)
  }
  
  private def handleTerminated(terminee: ActorRef) = {
    monitoredActors = monitoredActors - terminee
    loadMetrics = loadMetrics.filterKeys({
      case ActorAndMetric(actor, _) => actor != terminee
    })
  }
  
  private def highLoadMetricsForLogging = {
    loadMetrics.collect({
      case (ActorAndMetric(_, metricPath), HighLoad) => metricPath.head
    }).toSet.mkString(", ")
  }
}
