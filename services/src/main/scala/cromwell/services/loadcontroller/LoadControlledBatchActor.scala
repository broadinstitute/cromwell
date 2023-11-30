package cromwell.services.loadcontroller

import akka.actor.ActorRef
import cats.data.NonEmptyList
import cromwell.core.LoadConfig
import cromwell.core.actor.BatchActor
import cromwell.services.loadcontroller.LoadControlledBatchActor._
import cromwell.services.loadcontroller.LoadControllerService.{HighLoad, LoadMetric, NormalLoad}

object LoadControlledBatchActor {
  case object LoadControlledBatchActorTimerKey
  case object LoadControlledBatchActorTimerAction
}

/**
  * Can be mixed in a BatchActor to provide automated monitoring of the queue size based on a threshold value
  */
trait LoadControlledBatchActor[C] { this: BatchActor[C] =>
  def threshold: Int
  def serviceRegistryActor: ActorRef
  private val path =
    if (routed) NonEmptyList.of(context.parent.path.name, self.path.name) else NonEmptyList.one(self.path.name)

  timers.startSingleTimer(LoadControlledBatchActorTimerKey,
                          LoadControlledBatchActorTimerAction,
                          LoadConfig.MonitoringFrequency
  )

  private def weightToLoad(weight: Int) = if (weight > threshold) HighLoad else NormalLoad

  /**
    * Don't forget to chain this into your receive method to monitor the queue size:
    * override def receive = loadControlReceive.orElse(super.receive)
    */
  protected def loadControlReceive: Receive = { case LoadControlledBatchActorTimerAction =>
    serviceRegistryActor ! LoadMetric(path, weightToLoad(stateData.weight))
    timers.startSingleTimer(LoadControlledBatchActorTimerKey,
                            LoadControlledBatchActorTimerAction,
                            LoadConfig.MonitoringFrequency
    )
  }
}
