package cromwell.engine.workflow.tokens

import java.util.UUID
import java.util.concurrent.atomic.AtomicBoolean

import cromwell.core.JobExecutionToken
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.UnhoggableTokenPool._
import io.circe.generic.JsonCodec
import io.github.andrebeat.pool._

import scala.collection.immutable.HashSet
import scala.collection.mutable

final class UnhoggableTokenPool(val tokenType: JobExecutionTokenType) extends SimplePool[JobExecutionToken](
  capacity = tokenType.maxPoolSize.getOrElse(UnhoggableTokenPool.MaxCapacity),
  referenceType = ReferenceType.Strong,
  _factory = () => JobExecutionToken(tokenType, UUID.randomUUID()),
  _reset = Function.const(()),
  _dispose = Function.const(()),
  _healthCheck = Function.const(true)) {

  lazy val hogLimitOption: Option[Int] = tokenType match {
    case JobExecutionTokenType(_, Some(limit), hogFactor) if hogFactor > 1 =>
      Option(math.max(1, math.round(limit.floatValue() / hogFactor.floatValue())))
    case _ => None
  }

  private[this] val hogGroupAssignments: mutable.Map[String, HashSet[JobExecutionToken]] = mutable.Map.empty

  override def tryAcquire(): Option[Lease[JobExecutionToken]] = throw new UnsupportedOperationException("Use tryAcquire(hogGroup)")

  def available(hogGroup: String): UnhoggableTokenPoolAvailability = {

    hogLimitOption match {
      case None if leased < capacity => TokensAvailable
      case None => TokenTypeExhausted
      case Some(hogLimit) =>
        if (leased < capacity) {
          synchronized {
            if (hogGroupAssignments.get(hogGroup).forall(_.size < hogLimit)) {
              TokensAvailable
            } else {
              HogLimitExceeded
            }
          }
        } else TokenTypeExhausted
    }
  }

  def tryAcquire(hogGroup: String): UnhoggableTokenPoolResult = {

    hogLimitOption match {
      case Some(hogLimit) =>
        synchronized {
          val thisHogSet = hogGroupAssignments.getOrElse(hogGroup, HashSet.empty)

          if (thisHogSet.size < hogLimit) {
            super.tryAcquire() match {
              case Some(lease) =>
                val hoggingLease = new TokenHoggingLease(lease, hogGroup, this)
                hogGroupAssignments += hogGroup -> (thisHogSet + hoggingLease.get)
                hoggingLease
              case None => TokenTypeExhausted
            }
          } else {
            if (leased == capacity) TokenTypeExhausted else HogLimitExceeded
          }
        }
      case None =>
        super.tryAcquire() match {
          case Some(lease) =>
            val hoggingLease = new TokenHoggingLease(lease, hogGroup, this)
            hoggingLease
          case None => TokenTypeExhausted
        }
    }
  }

  def unhog(hogGroup: String, lease: Lease[JobExecutionToken]): Unit = {
    hogLimitOption foreach { _ =>
      synchronized {
        val newAssignment = hogGroupAssignments.getOrElse(hogGroup, HashSet.empty) - lease.get

        if (newAssignment.isEmpty) {
          hogGroupAssignments -= hogGroup
        } else {
          hogGroupAssignments += hogGroup -> newAssignment
        }
      }
    }
  }

  def poolState: TokenPoolState = {
    val (hogGroupUsages, hogLimitValue): (Option[Set[HogGroupState]], Option[Int]) = hogLimitOption match {
      case Some(hogLimit) =>
        synchronized {
          val entries: Set[HogGroupState] = hogGroupAssignments.toSet[(String, HashSet[JobExecutionToken])].map { case (hogGroup, set) =>
            HogGroupState(hogGroup, set.size, !hogGroupAssignments.get(hogGroup).forall(_.size < hogLimit))
          }
          (Option(entries), Option(hogLimit))
        }
      case None => (None, None)
    }

    TokenPoolState(hogGroupUsages, hogLimitValue, capacity, leased, leased < capacity)
  }
}

object UnhoggableTokenPool {
  // 10 Million
  val MaxCapacity = 10 * 1000 * 1000

  sealed trait UnhoggableTokenPoolResult

  final class TokenHoggingLease(lease: Lease[JobExecutionToken], hogGroup: String, pool: UnhoggableTokenPool) extends Lease[JobExecutionToken] with UnhoggableTokenPoolResult {
    private[this] val dirty = new AtomicBoolean(false)
    override protected[this] def a: JobExecutionToken = lease.get

    override protected[this] def handleRelease(): Unit = {
      if (dirty.compareAndSet(false, true)) {
        pool.unhog(hogGroup, lease)
      }
      lease.release()
    }
    override protected[this] def handleInvalidate(): Unit = {
      handleRelease()
      lease.invalidate()
    }
  }

  sealed trait UnhoggableTokenPoolAvailability { def available: Boolean }

  case object TokensAvailable extends UnhoggableTokenPoolAvailability {
    override def available: Boolean = true
  }

  /**
    * You didn't get a lease because the pool is empty
    */
  case object TokenTypeExhausted extends UnhoggableTokenPoolAvailability with UnhoggableTokenPoolResult {
    override def available: Boolean = false
  }

  /**
    * You didn't get a lease... because you're being a hog
    */
  case object HogLimitExceeded extends UnhoggableTokenPoolAvailability with UnhoggableTokenPoolResult {
    override def available: Boolean = false
  }

  @JsonCodec
  final case class HogGroupState(hogGroup: String, used: Int, atLimit: Boolean)

  @JsonCodec
  final case class TokenPoolState(hogGroups: Option[Set[HogGroupState]], hogLimit: Option[Int], capacity: Int, leased: Int, available: Boolean)

}
