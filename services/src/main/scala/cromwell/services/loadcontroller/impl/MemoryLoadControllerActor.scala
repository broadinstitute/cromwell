package cromwell.services.loadcontroller.impl

import akka.actor.{Actor, ActorRef, Props, Timers}

import scala.concurrent.duration._
import MemoryLoadControllerActor._
import cats.data.NonEmptyList
import cromwell.services.instrumentation.CromwellInstrumentationActor
import cromwell.services.loadcontroller.LoadControllerService.{HighLoad, LoadLevel, LoadMetric, NormalLoad}
object MemoryLoadControllerActor {
  case object MemoryLoadControlTimerKey
  case object MemoryLoadControlTimerAction
  private val runtime = Runtime.getRuntime
  private val FreeMemoryInstrumentationPath = NonEmptyList.one("freeMemoryInMB")
  private val MB = 1024L * 1024L
  private def freeMemoryInMB = (runtime.maxMemory() - (runtime.totalMemory() - runtime.freeMemory())) / MB
  def props(threshold: Long, frequency: FiniteDuration, serviceActor: ActorRef) = {
    Props(new MemoryLoadControllerActor(threshold, frequency, serviceActor))
  }
}

class MemoryLoadControllerActor(thresholdInMB: Long,
                                frequency: FiniteDuration,
                                override val serviceRegistryActor: ActorRef) extends Actor with Timers with CromwellInstrumentationActor {
  private var amortizedLoad: LoadLevel = NormalLoad
  private val nbRecordings = 10
  private var index: Int = 0
  private val loadRecordings = Array.fill[LoadLevel](nbRecordings) { NormalLoad }

  override def preStart() = {
    timers.startSingleTimer(MemoryLoadControlTimerKey, MemoryLoadControlTimerAction, frequency)
    super.preStart()
  }

  override def receive = {
    case MemoryLoadControlTimerAction => updateAndCheckMemory()
  }

  private [impl] def getFreeMemory = freeMemoryInMB

  private def updateAndCheckMemory() = {
    val currentFreeMemory = getFreeMemory
    val currentLoad = if (currentFreeMemory < thresholdInMB) HighLoad else NormalLoad
    // Update
    loadRecordings.update(index, currentLoad)
    index = (index + 1) % nbRecordings

    // Check if the load has been homogeneous for the past 10 checks
    loadRecordings.toSet.size match {
      // If yes and it's different from the current amortized load, update it
      case 1 if loadRecordings.head != amortizedLoad => amortizedLoad = loadRecordings.head
      case _ =>
    }

    serviceRegistryActor ! LoadMetric("Memory", amortizedLoad)
    sendGauge(FreeMemoryInstrumentationPath, currentFreeMemory)
    timers.startSingleTimer(MemoryLoadControlTimerKey, MemoryLoadControlTimerAction, frequency)
  }
}
