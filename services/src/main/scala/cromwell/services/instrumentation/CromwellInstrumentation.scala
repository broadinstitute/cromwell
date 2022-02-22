package cromwell.services.instrumentation

import java.time.{OffsetDateTime, Duration => JDuration}
import java.util.concurrent.TimeUnit

import akka.actor.{Actor, ActorRef, Timers}
import akka.dispatch.ControlMessage
import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.services.instrumentation.CromwellInstrumentation._
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.implicitConversions

object CromwellInstrumentation {

  val InstrumentationRate = ConfigFactory.load()
    .getConfig("system")
    .as[Option[FiniteDuration]]("instrumentation-rate")
    .getOrElse(5.seconds)

  implicit class InstrumentationPath (val internalPath: NonEmptyList[Either[String, (String, String)]]) {
    def :+(part: String): InstrumentationPath = internalPath.append(Left(part))
    def withParts(parts: String *): InstrumentationPath = withParts(parts.toList)
    def withParts(parts: List[String]): InstrumentationPath = internalPath.concat(parts.map(Left(_)))
    def withHighVariantPart(label: String, part: String): InstrumentationPath = withHighVariantPart(label -> part)
    def withHighVariantPart(tuple: (String, String)): InstrumentationPath = internalPath.append(Right(tuple))
    def withStatusCodeFailure(code: Option[Int]): InstrumentationPath = code match {
      case Some(value) => withHighVariantPart("code", value.toString)
      case None => this
    }
    def withThrowable(failure: Throwable, statusCodeExtractor: Throwable => Option[Int]): InstrumentationPath =
      internalPath.withStatusCodeFailure(statusCodeExtractor(failure))

    def concat(other: InstrumentationPath): InstrumentationPath = internalPath.concatNel(other.internalPath)

    /**
     * Get all path parts, in the order they were added, without special handling for high variant parts.
     * The "labels" recorded for high variant parts are unused.
     * @return a NeL of the parts of the instrumentation path
     */
    def getPath: NonEmptyList[String] = internalPath.map {
      case Left(part) => part
      case Right((_, part)) => part
    }

    /**
     * A best-effort method to get a non-empty list of path parts with the fewest number of high variant parts
     * possible. Any high variant parts not included in the list are returned as a label-part map.
     * @return a NeL of the instrumentation path, with as many high variant parts as possible instead returned in a map
     */
    def getPathLowVariants: (NonEmptyList[String], Map[String, String]) = {
      var nameParts = internalPath.collect { case Left(p) => p }
      var labelParts = internalPath.collect { case Right(p) => p }
      // path is a NeL, so if nameParts is empty then labelParts is not
      if (nameParts.isEmpty) {
        nameParts = labelParts.take(1).map(_._2)
        labelParts = labelParts.drop(1)
      }
      (NonEmptyList.fromListUnsafe(nameParts), labelParts.toMap)
    }
  }

  object InstrumentationPath {
    def withParts(part: String, additional: String *): InstrumentationPath = NonEmptyList.of(Left(part), additional.map(Left(_)):_*)
    def withHighVariantPart(tuple: (String, String)): InstrumentationPath = NonEmptyList.of(Right(tuple))
    def withHighVariantPart(label: String, part: String): InstrumentationPath = withHighVariantPart(label -> part)
  }

  @deprecated("CromwellInstrumentation no longer exposes NeLs; this method will be removed", "Cromwell 75")
  implicit def stringToNel(str: String): NonEmptyList[String] = NonEmptyList.of(str)
}

trait CromwellInstrumentationActor extends CromwellInstrumentation { this: Actor =>
  override def instrumentationSender: ActorRef = self
}

trait CromwellInstrumentation {
  protected def instrumentationSender: ActorRef = ActorRef.noSender

  def serviceRegistryActor: ActorRef

  /**
    * Builds a bucket from path.
    * The cromwell bucket prefix is always prepended:
    * cromwell.[prefix].path
    */
  final private def makeBucket(path: InstrumentationPath, prefix: Option[String]): CromwellBucket = {
    CromwellBucket(prefix.toList, path)
  }

  /**
    * Creates an increment message for the given bucket
    */
  private final def countMessage(path: InstrumentationPath, count: Long, prefix: Option[String]): InstrumentationServiceMessage = {
    InstrumentationServiceMessage(CromwellCount(makeBucket(path, prefix), count))
  }

  /**
    * Increment the counter for the given bucket
    */
  protected final def count(path: InstrumentationPath, count: Long, prefix: Option[String] = None): Unit = {
    serviceRegistryActor.tell(countMessage(path, count, prefix), instrumentationSender)
  }

  /**
    * Creates an increment message for the given bucket
    */
  private final def incrementMessage(path: InstrumentationPath, prefix: Option[String]): InstrumentationServiceMessage = {
    InstrumentationServiceMessage(CromwellIncrement(makeBucket(path, prefix)))
  }

  /**
    * Increment the counter for the given bucket
    */
  protected final def increment(path: InstrumentationPath, prefix: Option[String] = None): Unit = {
    serviceRegistryActor.tell(incrementMessage(path, prefix), instrumentationSender)
  }

  /**
    * Creates a gauge message for the given bucket
    */
  private final def gaugeMessage(path: InstrumentationPath, value: Long, prefix: Option[String]) = {
    InstrumentationServiceMessage(CromwellGauge(makeBucket(path, prefix), value))
  }

  /**
    * Set the bucket to the gauge value
    */
  protected final def sendGauge(path: InstrumentationPath, value: Long, prefix: Option[String] = None): Unit = {
    serviceRegistryActor.tell(gaugeMessage(path, value, prefix), instrumentationSender)
  }

  /**
    * Creates a timing message for the given bucket and duration
    */
  private final def timingMessage(path: InstrumentationPath, duration: FiniteDuration, prefix: Option[String]) = {
    InstrumentationServiceMessage(CromwellTiming(makeBucket(path, prefix), duration))
  }

  /**
    * Add a timing information for the given bucket
    */
  protected final def sendTiming(path: InstrumentationPath, duration: FiniteDuration, prefix: Option[String] = None) = {
    serviceRegistryActor.tell(timingMessage(path, duration, prefix), instrumentationSender)
  }

  def calculateTimeDifference(startTime: OffsetDateTime, endTime: OffsetDateTime): FiniteDuration = {
    FiniteDuration(JDuration.between(startTime, endTime).toMillis, TimeUnit.MILLISECONDS)
  }
  def calculateTimeSince(startTime: OffsetDateTime): FiniteDuration = calculateTimeDifference(startTime, OffsetDateTime.now())
}

/**
  * Helper trait to provide a scheduler function that can be used for instrumentation purposes
  */
trait CromwellInstrumentationScheduler { this: Actor with Timers =>
  private case object InstrumentationTimerKey
  private case object InstrumentationTimerAction extends ControlMessage

  def startInstrumentationTimer() = {
    timers.startSingleTimer(InstrumentationTimerKey, InstrumentationTimerAction, InstrumentationRate)
  }

  protected def instrumentationReceive(instrumentationAction: () => Unit): Receive = {
    case InstrumentationTimerAction =>
      instrumentationAction()
      timers.startSingleTimer(InstrumentationTimerKey, InstrumentationTimerAction, InstrumentationRate)
  }
}
