package cromwell.engine.workflow.tokens

import java.util.UUID
import java.util.concurrent.atomic.AtomicBoolean

import cromwell.core.JobExecutionToken
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.UnhoggableTokenPool._
import io.github.andrebeat.pool._
import spray.json._

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
    case JobExecutionTokenType(_, None, _) => None
    case JobExecutionTokenType(_, Some(limit), hogFactor) if hogFactor > 1 => Option(math.max(1, math.round(limit.floatValue() / hogFactor.floatValue())))
    case JobExecutionTokenType(_, _, _) => None
  }

  private[this] val hogGroupAssignments: mutable.Map[String, HashSet[JobExecutionToken]] = mutable.Map.empty

  override def tryAcquire(): Option[Lease[JobExecutionToken]] = throw new UnsupportedOperationException("Use tryAcquire(hogGroup)")

  def available(hogGroup: String): UnhoggableTokenPoolAvailability = {

    hogLimitOption match {
      case None if leased < capacity => TokensAvailable
      case None => ComeBackLater
      case Some(hogLimit) =>
        if (leased < capacity) {
          synchronized {
            if (hogGroupAssignments.get(hogGroup).forall(_.size < hogLimit)) {
              TokensAvailable
            } else {
              Oink
            }
          }
        } else ComeBackLater
    }
  }

  def tryAcquire(hogGroup: String): UnhoggableTokenPoolResult = {

    hogLimitOption match {
      case Some(hogLimit) =>
        synchronized {
          val thisHogSet = hogGroupAssignments.getOrElse(hogGroup, HashSet.empty)

          if (thisHogSet.size + 1 <= hogLimit) {
            super.tryAcquire() match {
              case Some(lease) =>
                val hoggingLease = new TokenHoggingLease(lease, hogGroup, this)
                hogGroupAssignments += hogGroup -> (thisHogSet + hoggingLease.get)
                hoggingLease
              case None => ComeBackLater
            }
          } else {
            if (leased == capacity) ComeBackLater else Oink
          }
        }
      case None =>
        super.tryAcquire() match {
          case Some(lease) =>
            val hoggingLease = new TokenHoggingLease(lease, hogGroup, this)
            hoggingLease
          case None => ComeBackLater
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

  def poolState: JsObject = {
    val (hogGroupUsages, hogLimitValue): (JsValue, JsValue) = hogLimitOption match {
      case Some(hogLimit) =>
        val entries = hogGroupAssignments.map { case (hogGroup, set) =>
          JsObject(Map(
            "hog group" -> JsString(hogGroup),
            "used" -> JsNumber(set.size),
            "available" -> JsBoolean(hogGroupAssignments.get(hogGroup).forall(_.size < hogLimit))
          ))
        }
        (JsArray(entries.toVector), JsNumber(hogLimit))
      case None => (JsNull, JsNull)
    }

    JsObject(Map(
      "hog groups" -> hogGroupUsages,
      "hog limit" -> hogLimitValue,
      "capacity" -> JsNumber(capacity),
      "leased" -> JsNumber(leased)
    ))
  }
}

object UnhoggableTokenPool {
  // 10 Million
  val MaxCapacity = 10 * 1000 * 1000

  trait UnhoggableTokenPoolResult

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

  trait UnhoggableTokenPoolAvailability { def available: Boolean }

  case object TokensAvailable extends UnhoggableTokenPoolAvailability {
    override def available: Boolean = true
  }

  /**
    * You didn't get a lease because the pool is empty
    */
  case object ComeBackLater extends UnhoggableTokenPoolAvailability with UnhoggableTokenPoolResult {
    override def available: Boolean = false
  }

  /**
    * You didn't get a lease... because you're being a hog
    */
  case object Oink extends UnhoggableTokenPoolAvailability with UnhoggableTokenPoolResult {
    override def available: Boolean = false
  }

}
