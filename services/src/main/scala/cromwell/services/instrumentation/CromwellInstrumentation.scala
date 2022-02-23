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

object CromwellInstrumentation {

  val InstrumentationRate = ConfigFactory.load()
    .getConfig("system")
    .as[Option[FiniteDuration]]("instrumentation-rate")
    .getOrElse(5.seconds)

  /**
   * A representation of a metric name, with different output formats for different kinds of instrumentation services.
   * Those designed to handle one time series per metric name can get a simple list of name parts with [[getFlatPath]],
   * while those designed to handle multiple time series per metric name can get a shorter metric name and
   * frequently-changing, high-variant "labels" separately via [[getPathAndLabels]].
   *
   * This distinction exists because instrumentation services meant to handle multiple time series per metric name
   * (like Prometheus) are poorly built to store numerous metric names or query across them. The more complex
   * [[getPathAndLabels]] output allows the backends for those instrumentation services to "play nice" with their
   * associated ecosystems.
   *
   * @param internalPath the runtime representation of this class, where the 'left' is a simple, always-present
   *                     low-variant name part, and the 'right' is a key-value pair for a frequently-changing
   *                     high-variant name part
   */
  implicit class InstrumentationPath (val internalPath: NonEmptyList[Either[String, (String, String)]]) extends AnyVal {
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
      withStatusCodeFailure(statusCodeExtractor(failure))

    def concat(other: InstrumentationPath): InstrumentationPath = internalPath.concatNel(other.internalPath)

    /**
     * Get all path parts as an ordered list, handling high-variant parts by extracting only the value of the key-value
     * pair.
     * @return a NeL of the parts of the instrumentation path
     */
    def getFlatPath: NonEmptyList[String] = internalPath.map {
      case Left(part) => part
      case Right((_, part)) => part
    }

    /**
     * Get path parts as an ordered list, with as many high-variant parts as possible excluded and instead returned in a
     * label map (the invariant is that the path part list not be empty, so this method will handle the first
     * high-variant part like [[getFlatPath]] if there are zero normal parts).
     * @return a NeL of the parts of the instrumentation path, and a separate map for labels
     */
    def getPathAndLabels: (NonEmptyList[String], Map[String, String]) = {
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

  /**
   * Companion object for [[InstrumentationPath]], containing constructor methods to help maintain the non-empty path
   * invariant.
   */
  object InstrumentationPath {
    def withParts(part: String, additional: String *): InstrumentationPath = NonEmptyList.of(Left(part), additional.map(Left(_)):_*)
    def withHighVariantPart(tuple: (String, String)): InstrumentationPath = NonEmptyList.of(Right(tuple))
    def withHighVariantPart(label: String, part: String): InstrumentationPath = withHighVariantPart(label -> part)
  }
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
