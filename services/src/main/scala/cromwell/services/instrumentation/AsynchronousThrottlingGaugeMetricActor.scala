package cromwell.services.instrumentation

import akka.actor.{ActorRef, LoggingFSM, Props, Status}
import akka.pattern.pipe
import cats.data.NonEmptyList
import cromwell.services.instrumentation.AsynchronousThrottlingGaugeMetricActor._

import scala.concurrent.{ExecutionContext, Future}

class AsynchronousThrottlingGaugeMetricActor(metricPath: NonEmptyList[String],
                                             instrumentationPrefix: Option[String],
                                             override val serviceRegistryActor: ActorRef)
  extends LoggingFSM[AsynchronousThrottlingGaugeMetricActorState, Unit]
    with CromwellInstrumentation {

  implicit val ec = context.dispatcher

  startWith(WaitingForMetricCalculationRequestOrMetricValue, ())

  when(WaitingForMetricCalculationRequestOrMetricValue) {
    case Event(MetricValue(value), _) =>
      sendGaugeAndStay(value)
    case Event(CalculateMetricValue(calculateMetricValueFunction), _) =>
      calculateMetricValueFunction(ec)
        .map(intValue => FinishedMetricValueCalculation(intValue.toLong))
        .pipeTo(self)
      goto(MetricCalculationInProgress)
  }

  when(MetricCalculationInProgress) {
    case Event(CalculateMetricValue(_), _) =>
      // metric actor is already busy calculating metric value, so we dismiss this request
      stay()
    case Event(MetricValue(value), _) =>
      // pass through the pre-calculated metric value and stay, continuing to wait for ongoing metric value calculation to finish
      sendGaugeAndStay(value)
    case Event(FinishedMetricValueCalculation(value), _) =>
      sendGauge(metricPath, value, instrumentationPrefix)
      goto(WaitingForMetricCalculationRequestOrMetricValue)
    case Event(Status.Failure(ex), _) =>
      log.error(s"Cannot send gauge metric value: error occurred while evaluating Future value.", ex)
      goto(WaitingForMetricCalculationRequestOrMetricValue)
  }

  whenUnhandled {
    case Event(unexpected, _) =>
      log.warning(s"Programmer error: this actor should not receive message $unexpected from ${sender.path} while in state $stateName")
      stay()
  }

  private def sendGaugeAndStay(metricValue: Long): State = {
    sendGauge(metricPath, metricValue, instrumentationPrefix)
    stay()
  }

}

object AsynchronousThrottlingGaugeMetricActor {

  def props(metricPath: NonEmptyList[String], instrumentationPrefix: Option[String], serviceRegistryActor: ActorRef) =
    Props(new AsynchronousThrottlingGaugeMetricActor(metricPath, instrumentationPrefix, serviceRegistryActor))

  // states
  sealed trait AsynchronousThrottlingGaugeMetricActorState
  case object WaitingForMetricCalculationRequestOrMetricValue extends AsynchronousThrottlingGaugeMetricActorState
  case object MetricCalculationInProgress extends AsynchronousThrottlingGaugeMetricActorState

  // messages
  case class CalculateMetricValue(calculateMetricValueFunction: ExecutionContext => Future[Int])
  case class MetricValue(value: Long)
  case class FinishedMetricValueCalculation(calculatedValue: Long)

}
