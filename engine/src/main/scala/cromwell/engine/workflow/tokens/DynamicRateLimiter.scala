package cromwell.engine.workflow.tokens

import akka.actor.{Actor, ActorLogging, Timers}
import akka.dispatch.ControlMessage
import cromwell.engine.workflow.tokens.DynamicRateLimiter._
import cromwell.services.loadcontroller.LoadControllerService.{HighLoad, NormalLoad}

import scala.concurrent.duration._

/**
  * Simply listen to load alerts and freeze the token distribution if the load is high, restore it when load is back to normal.
  */
trait DynamicRateLimiter { this: Actor with Timers with ActorLogging =>
  protected def distributionRate: Rate

  protected def ratePreStart(): Unit = {
    startDistributionTimer()
    log.info("{} - Distribution rate: {}.", self.path.name, distributionRate)
  }

  protected def rateReceive: Receive = {
    case ResetAction if distributionRate.n != 0 => releaseTokens()
    case HighLoad => highLoad()
    case NormalLoad => backToNormal()
  }

  private def startDistributionTimer() = if (!distributionRate.isZero && !timers.isTimerActive(ResetKey)) {
    timers.startPeriodicTimer(ResetKey, ResetAction, distributionRate.per)
  }

  private def releaseTokens() = {
    self ! TokensAvailable(distributionRate.n)
  }

  // When load is high, freeze token distribution
  private def highLoad(doLogging: Boolean = true) = {
    // This guarantees that even if the ResetAction message has already been sent and is in the message queue it will be removed in time
    timers.cancel(ResetKey)
    log.warning("{} - High load alert. Freeze token distribution.", self.path.name)
  }

  // When back to normal, restart the token distribution timer if needed
  private def backToNormal() = {
    startDistributionTimer()
    log.info("{} - Load back to normal", self.path.name)
  }
}

object DynamicRateLimiter {
  private case object ResetKey
  private case object ResetAction extends ControlMessage

  case class TokensAvailable(n: Int) extends ControlMessage

  implicit class RateEnhancer(val n: Int) extends AnyVal {
    def per(duration: FiniteDuration) = Rate(n, duration)
  }

  case class Rate(n: Int, per: FiniteDuration) {
    def isZero = n == 0 || per == Duration.Zero
    override def toString = s"$n per ${per.toSeconds} seconds"
  }
}
