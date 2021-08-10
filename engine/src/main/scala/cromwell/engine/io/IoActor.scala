package cromwell.engine.io

import akka.NotUsed
import akka.actor.{Actor, ActorLogging, ActorRef, Props, Timers}
import akka.dispatch.ControlMessage
import akka.stream._
import akka.stream.scaladsl.{Flow, GraphDSL, Merge, Partition, Sink, Source, SourceQueueWithComplete}
import com.typesafe.config.Config
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.actor.StreamActorHelper
import cromwell.core.actor.StreamIntegration.StreamContext
import cromwell.core.io.Throttle._
import cromwell.core.io.{IoAck, IoCommand, Throttle}
import cromwell.core.{Dispatcher, LoadConfig}
import cromwell.engine.instrumentation.IoInstrumentation
import cromwell.engine.io.IoActor._
import cromwell.engine.io.gcs.GcsBatchFlow.GcsBatchFlowConfig
import cromwell.engine.io.gcs.{GcsBatchCommandContext, ParallelGcsBatchFlow}
import cromwell.engine.io.nio.NioFlow
import cromwell.engine.io.nio.NioFlow.NioFlowConfig
import cromwell.filesystems.gcs.batch.GcsBatchIoCommand
import cromwell.services.loadcontroller.LoadControllerService.{HighLoad, LoadMetric, NormalLoad}
import net.ceedubs.ficus.readers.ValueReader

import java.time.OffsetDateTime
import java.time.temporal.ChronoUnit
import scala.concurrent.ExecutionContext
import scala.concurrent.duration._


/**
  * Actor that performs IO operations asynchronously using akka streams
  * 
  * @param ioConfig IoActor configuration class
  * @param materializer actor materializer to run the stream
  * @param serviceRegistryActor actorRef for the serviceRegistryActor
  */
final class IoActor(ioConfig: IoConfig,
                    override val serviceRegistryActor: ActorRef,
                    applicationName: String)(implicit val materializer: ActorMaterializer)
  extends Actor with ActorLogging with StreamActorHelper[IoCommandContext[_]] with IoInstrumentation with Timers {
  implicit val ec: ExecutionContext = context.dispatcher

  // IntelliJ disapproves of mutable state in Actors, but this should be safe as long as access occurs only in
  // the `receive` method. Alternatively IntelliJ does suggest a `become` workaround we might try in the future.

  //noinspection ActorMutableStateInspection
  private var backpressureExpiration: Option[OffsetDateTime] = None

  /**
    * Method for instrumentation to be executed when a IoCommand failed and is being retried.
    * Can be passed to flows so they can invoke it when necessary.
    */
  private def onRetry(commandContext: IoCommandContext[_])(throwable: Throwable): Unit = {
    incrementIoRetry(commandContext.request, throwable)
  }

  override def preStart(): Unit = {
    // On start up, let the controller know that the load is normal
    serviceRegistryActor ! LoadMetric("IO", NormalLoad)
    super.preStart()
  }

  private [io] lazy val defaultFlow =
    new NioFlow(
      parallelism = ioConfig.nio.parallelism,
      onRetryCallback = onRetry,
      onBackpressure = onBackpressure,
      numberOfAttempts = ioConfig.numberOfAttempts,
      commandBackpressureStaleness = ioConfig.commandBackpressureStaleness)
      .flow
      .withAttributes(ActorAttributes.dispatcher(Dispatcher.IoDispatcher))

  private [io] lazy val gcsBatchFlow =
    new ParallelGcsBatchFlow(
      config = ioConfig.gcsBatch,
      scheduler = context.system.scheduler,
      onRetry = onRetry,
      onBackpressure = onBackpressure,
      applicationName = applicationName,
      commandBackpressureStaleness = ioConfig.commandBackpressureStaleness)
      .flow
      .withAttributes(ActorAttributes.dispatcher(Dispatcher.IoDispatcher))

  private val source = Source.queue[IoCommandContext[_]](ioConfig.queueSize, OverflowStrategy.dropNew)

  private val flow = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._

    val input = builder.add(Flow[IoCommandContext[_]])

    // Partitions requests between gcs batch, and single nio requests
    val batchPartitioner = builder.add(Partition[IoCommandContext[_]](2, {
      case _: GcsBatchCommandContext[_, _] => 0
      case _ => 1
    }))

    // Sub flow for batched gcs requests
    val batches = batchPartitioner.out(0) collect { case batch: GcsBatchCommandContext[_, _] => batch }

    // Sub flow for single nio requests
    val defaults = batchPartitioner.out(1) collect { case default: DefaultCommandContext[_] => default }

    // Merge results from both flows back together
    val merger = builder.add(Merge[IoResult](2))

    // Flow processing nio requests
    val defaultFlowPorts = builder.add(defaultFlow)

    // Flow processing gcs batch requests
    val batchFlowPorts = builder.add(gcsBatchFlow)

    input ~> batchPartitioner
             defaults.outlet ~> defaultFlowPorts ~> merger
             batches.outlet ~> batchFlowPorts ~> merger

    FlowShape[IoCommandContext[_], IoResult](input.in, merger.out)
  }

  private val throttledFlow = ioConfig.throttle map { t =>
    Flow[IoCommandContext[_]]
      .throttle(t.elements, t.per, t.maximumBurst, ThrottleMode.Shaping)
      .via(flow)
  } getOrElse flow

  private val instrumentationSink = Sink.foreach[IoResult](instrumentIoResult)

  override protected lazy val streamSource: Source[
    (Any, IoCommandContext[_]),
    SourceQueueWithComplete[IoCommandContext[_]]
  ] = source
    .via(throttledFlow)
    .alsoTo(instrumentationSink)
    .withAttributes(ActorAttributes.dispatcher(Dispatcher.IoDispatcher))

  override def onBackpressure(scale: Option[Double] = None): Unit = {
    incrementBackpressure()
    serviceRegistryActor ! LoadMetric("IO", HighLoad)

    val uncappedDelay = scale.getOrElse(1.0d) * LoadConfig.IoNormalWindowMinimum
    val cappedDelay = FiniteDuration(LoadConfig.IoNormalWindowMaximum.min(uncappedDelay).toMillis, MILLISECONDS)

    self ! BackPressure(cappedDelay)
  }

  override def actorReceive: Receive = {
    /* GCS Batch command with context */
    case (clientContext: Any, gcsBatchCommand: GcsBatchIoCommand[_, _]) =>
      val replyTo = sender()
      val commandContext = GcsBatchCommandContext(
        request = gcsBatchCommand,
        maxAttemptsNumber = ioConfig.numberOfAttempts,
        replyTo = replyTo,
        clientContext = Option(clientContext))
      sendToStream(commandContext)

    /* GCS Batch command without context */
    case gcsBatchCommand: GcsBatchIoCommand[_, _] =>
      val replyTo = sender()
      val commandContext = GcsBatchCommandContext(
        request = gcsBatchCommand,
        maxAttemptsNumber = ioConfig.numberOfAttempts,
        replyTo = replyTo
      )
      sendToStream(commandContext)

    /* Default command with context */
    case (clientContext: Any, command: IoCommand[_]) =>
      val replyTo = sender()
      val commandContext = DefaultCommandContext(command, replyTo, Option(clientContext))
      sendToStream(commandContext)

    /* Default command without context */
    case command: IoCommand[_] =>
      val replyTo = sender()
      val commandContext = DefaultCommandContext(command, replyTo)
      sendToStream(commandContext)

    case BackPressure(duration) =>
      lazy val proposedExpiry = OffsetDateTime.now().plusNanos(duration.toNanos)

      // Because this method will be called every time we backpressure, the timer will be overridden every
      // time until we're not backpressuring anymore
      backpressureExpiration match {
        case None =>
          // Start a new backpressure
          val durationSeconds = duration.toMillis / 1000.0
          log.info("Beginning IoActor backpressure for {} seconds", f"$durationSeconds%.2f")
          timers.startSingleTimer(BackPressureTimerResetKey, BackPressureTimerResetAction, duration)
          backpressureExpiration = Option(proposedExpiry)

        case Some(expiry) if expiry.isBefore(proposedExpiry) =>
          // Extend the current backpressure
          val extension = expiry.until(proposedExpiry, ChronoUnit.MILLIS)
          // There can be a lot of very short extensions when bursts of I/O requests are rejected from a full I/O
          // queue. Don't clutter the logs with messages for these.
          if (extension > ioConfig.backPressureExtensionLogThreshold.toMillis) {
            val extensionSeconds = extension / 1000.0
            log.info("Extending IoActor backpressure {} seconds", f"$extensionSeconds%.2f")
          }

          val newExpiration = OffsetDateTime.now().until(proposedExpiry, ChronoUnit.MILLIS)
          timers.startSingleTimer(BackPressureTimerResetKey, BackPressureTimerResetAction, FiniteDuration(newExpiration, MILLISECONDS))
          backpressureExpiration = Option(proposedExpiry)

        case _ => // Ignore proposed expiries that would be before the current expiry
      }

    case BackPressureTimerResetAction =>
      log.info("IoActor backpressure off")
      backpressureExpiration = None
      serviceRegistryActor ! LoadMetric("IO", NormalLoad)
  }
}

trait IoCommandContext[T] extends StreamContext {
  val creationTime: OffsetDateTime = OffsetDateTime.now()
  def request: IoCommand[T]
  def replyTo: ActorRef
  def fail(failure: Throwable): IoResult = (request.fail(failure), this)
  def failReadForbidden(failure: Throwable, forbiddenPath: String): IoResult = (request.failReadForbidden(failure, forbiddenPath), this)
  def success(value: T): IoResult = (request.success(value), this)
}

object IoActor {

  import net.ceedubs.ficus.Ficus._

  /** Flow that can consume an IoCommandContext and produce an IoResult */
  type IoFlow = Flow[IoCommandContext[_], IoResult, NotUsed]

  /** Result type of an IoFlow, contains the original command context and the final IoAck response. */
  type IoResult = (IoAck[_], IoCommandContext[_])

  case class DefaultCommandContext[T](request: IoCommand[T], replyTo: ActorRef, override val clientContext: Option[Any] = None) extends IoCommandContext[T]

  case object BackPressureTimerResetKey

  case object BackPressureTimerResetAction extends ControlMessage

  case class BackPressure(duration: FiniteDuration) extends ControlMessage

  def props(ioConfig: IoConfig,
            serviceRegistryActor: ActorRef,
            applicationName: String,
           )
           (implicit materializer: ActorMaterializer): Props = {
    Props(new IoActor(ioConfig, serviceRegistryActor, applicationName)).withDispatcher(IoDispatcher)
  }

  case class IoConfig(queueSize: Int,
                      numberOfAttempts: Int,
                      commandBackpressureStaleness: FiniteDuration,
                      backPressureExtensionLogThreshold: FiniteDuration,
                      ioNormalWindowMinimum: FiniteDuration,
                      ioNormalWindowMaximum: FiniteDuration,
                      nio: NioFlowConfig,
                      gcsBatch: GcsBatchFlowConfig,
                      throttle: Option[Throttle])

  implicit val ioConfigReader: ValueReader[IoConfig] = (config: Config, _: String) => {

    val loadControl: Config = config.as[Config]("load-control")
    val queueSize: Int = loadControl.as[Int]("io-queue-size")
    val ioNormalWindowMinimum: FiniteDuration = loadControl.as[FiniteDuration]("io-normal-window-minimum")
    val ioNormalWindowMaximum: FiniteDuration = loadControl.as[FiniteDuration]("io-normal-window-maximum")

    val io: Config = config.as[Config]("system.io")
    val nioConfig = io.as[NioFlowConfig]("nio")
    val gcsConfig = io.as[GcsBatchFlowConfig]("gcs")
    val commandBackpressureStaleness = io.as[FiniteDuration]("command-backpressure-staleness")
    val backpressureExtensionLogThreshold = io.as[FiniteDuration]("backpressure-extension-log-threshold")
    val numberOfAttempts: Int = io.as[Int]("number-of-attempts")

    val throttle = io.as[Option[Throttle]]("throttle")

    IoConfig(
      queueSize = queueSize,
      numberOfAttempts = numberOfAttempts,
      commandBackpressureStaleness = commandBackpressureStaleness,
      backPressureExtensionLogThreshold = backpressureExtensionLogThreshold,
      ioNormalWindowMinimum = ioNormalWindowMinimum,
      ioNormalWindowMaximum = ioNormalWindowMaximum,
      throttle = throttle,
      nio = nioConfig,
      gcsBatch = gcsConfig
    )
  }
}
