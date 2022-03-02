package cromwell.services.instrumentation

import java.time.{OffsetDateTime, Duration => JDuration}
import java.util.concurrent.TimeUnit
import akka.actor.{Actor, ActorRef, Timers}
import akka.dispatch.ControlMessage
import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.services.instrumentation.CromwellInstrumentation.InstrumentationPath.requireNotEmpty
import cromwell.services.instrumentation.CromwellInstrumentation._
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._

object CromwellInstrumentation {

  val InstrumentationRate: FiniteDuration = ConfigFactory.load()
    .getConfig("system")
    .as[Option[FiniteDuration]]("instrumentation-rate")
    .getOrElse(5.seconds)

  /**
   * A part of a metric name that does not take on many other values. These low variant parts will always be included
   * in a final metric name regardless of the backing instrumentation service.
   * An empty [[LowVariantPart]] should be entirely excluded from metric names.
   */
  private type LowVariantPart = String

  /**
   * A part of a metric name that may take on many other values, represented by a key-value pair. The key should be
   * consistent for a given type of part (like `"code" -> response.status.intValue.toString`) to facilitate querying in
   * some backing instrumentation services. These high variant parts will either be included in the metric name directly
   * of will be supplied alongside as labels.
   * A [[HighVariantPart]] with an empty value should be entirely excluded from any metric names but should be included
   * in any labels to keep the label keys consistent. Empty keys are not permitted.
   */
  private type HighVariantPart = (String, String)

  /**
   * The lossless internal representation of a metric name, maintaining both part ordering and part variance.
   * The leading part should not be an empty string so as to conceptually maintain the non-empty invariant.
   */
  private type InternalPath = NonEmptyList[Either[LowVariantPart, HighVariantPart]]

  /**
   * A wrapping type around an [[InternalPath]] providing methods for building metric names and retrieving them in
   * formats suitable for different backing instrumentation services.
   *
   * @param internalPath the runtime representation of this class
   */
  implicit class InstrumentationPath (val internalPath: InternalPath) extends AnyVal {
    def :+(part: String): InstrumentationPath = internalPath.append(Left(part))
    def withParts(parts: String *): InstrumentationPath = internalPath.concat(parts.toList.map(Left(_)))
    def withHighVariantPart(label: String, part: String): InstrumentationPath = {
      requireNotEmpty("label", label)
      internalPath.append(Right(label -> part))
    }
    def withHighVariantPart(label: String, part: Option[String]): InstrumentationPath =
      withHighVariantPart(label, part.getOrElse(""))
    def withStatusCodeFailure(code: Option[Int]): InstrumentationPath =
      withHighVariantPart("code", code.map(_.toString))
    def withThrowable(failure: Throwable, statusCodeExtractor: Throwable => Option[Int]): InstrumentationPath =
      withStatusCodeFailure(statusCodeExtractor(failure))

    def concat(other: InstrumentationPath): InstrumentationPath = internalPath.concatNel(other.internalPath)

    /**
     * Get all path parts as an ordered list, handling [[HighVariantPart]] by extracting only the value of the key-value
     * pair.
     * @return a NeL of the parts of the instrumentation path
     */
    def getFlatPath: NonEmptyList[String] = NonEmptyList.fromListUnsafe(
      internalPath.map {
        case Left(part) => part
        case Right((_, part)) => part
      }
        // the leading part is never blank, per InstrumentationPath's object
        .filterNot(_.isBlank)
    )

    /**
     * Get path parts as an ordered list, with as many [[HighVariantPart]] as possible excluded and instead returned in
     * a label map (the invariant is that the path part list not be empty, so if there are no [[LowVariantPart]] then
     * the first [[HighVariantPart]] will be handled like [[getFlatPath]]).
     * @return a NeL of the parts of the instrumentation path, and a separate map for labels
     */
    def getPathAndLabels: (NonEmptyList[String], Map[String, String]) = {
      var nameParts = internalPath.collect { case Left(p) => p }
      var labelParts = internalPath.collect { case Right(p) => p }
      // internalPath is a NeL, so if nameParts is empty then labelParts is not
      if (nameParts.isEmpty) {
        nameParts = labelParts.take(1).map(_._2)
        labelParts = labelParts.drop(1)
      }
      // the leading part is never blank, per InstrumentationPath's object
      (NonEmptyList.fromListUnsafe(nameParts.filterNot(_.isBlank)), labelParts.toMap)
    }
  }

  /**
   * Companion object for [[InstrumentationPath]], containing constructor methods to help maintain the non-empty path
   * invariant (including that the leading part of a path cannot be empty).
   */
  object InstrumentationPath {
    private def requireNotEmpty(descriptor: String, value: String): Unit =
      require(!value.isBlank, s"InstrumentationPath $descriptor cannot be empty here")
    def withParts(part: String, additional: String *): InstrumentationPath = {
      requireNotEmpty("part", part)
      NonEmptyList.of(Left(part), additional.map(Left(_)):_*)
    }
    def withHighVariantPart(label: String, part: String): InstrumentationPath = {
      requireNotEmpty("label", label)
      requireNotEmpty("part", part)
      NonEmptyList.of(Right(label -> part))
    }
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
  protected final def sendTiming(path: InstrumentationPath, duration: FiniteDuration, prefix: Option[String] = None): Unit = {
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

  def startInstrumentationTimer(): Unit = {
    timers.startSingleTimer(InstrumentationTimerKey, InstrumentationTimerAction, InstrumentationRate)
  }

  protected def instrumentationReceive(instrumentationAction: () => Unit): Receive = {
    case InstrumentationTimerAction =>
      instrumentationAction()
      timers.startSingleTimer(InstrumentationTimerKey, InstrumentationTimerAction, InstrumentationRate)
  }
}
