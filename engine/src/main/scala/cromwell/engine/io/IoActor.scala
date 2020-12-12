package cromwell.engine.io

import java.time.OffsetDateTime

import akka.NotUsed
import akka.actor.{Actor, ActorLogging, ActorRef, Props, Timers}
import akka.dispatch.ControlMessage
import akka.stream._
import akka.stream.scaladsl.{Flow, GraphDSL, Merge, Partition, Sink, Source, SourceQueueWithComplete}
import com.typesafe.config.ConfigFactory
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.actor.StreamActorHelper
import cromwell.core.actor.StreamIntegration.StreamContext
import cromwell.core.io.{IoAck, IoCommand, Throttle}
import cromwell.core.{Dispatcher, LoadConfig}
import cromwell.engine.instrumentation.IoInstrumentation
import cromwell.engine.io.IoActor._
import cromwell.engine.io.gcs.{GcsBatchCommandContext, ParallelGcsBatchFlow}
import cromwell.engine.io.nio.NioFlow
import cromwell.filesystems.gcs.batch.GcsBatchIoCommand
import cromwell.services.loadcontroller.LoadControllerService.{HighLoad, LoadMetric, NormalLoad}

import scala.concurrent.ExecutionContext

/**
  * Actor that performs IO operations asynchronously using akka streams
  * 
  * @param queueSize size of the queue
  * @param throttle optional throttler to control the throughput of requests.
  *                 Applied to ALL incoming requests
  * @param materializer actor materializer to run the stream
  * @param serviceRegistryActor actorRef for the serviceRegistryActor
  */
final class IoActor(queueSize: Int,
                    nioParallelism: Int,
                    gcsParallelism: Int,
                    throttle: Option[Throttle],
                    override val serviceRegistryActor: ActorRef,
                    applicationName: String)(implicit val materializer: ActorMaterializer)
  extends Actor with ActorLogging with StreamActorHelper[IoCommandContext[_]] with IoInstrumentation with Timers {
  implicit val ec: ExecutionContext = context.dispatcher

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
    new NioFlow(parallelism = nioParallelism, onRetry)
      .flow
      .withAttributes(ActorAttributes.dispatcher(Dispatcher.IoDispatcher))
  private [io] lazy val gcsBatchFlow = new ParallelGcsBatchFlow(parallelism = gcsParallelism, batchSize = 100, context.system.scheduler, onRetry, applicationName).flow.withAttributes(ActorAttributes.dispatcher(Dispatcher.IoDispatcher))
  
  private val source = Source.queue[IoCommandContext[_]](queueSize, OverflowStrategy.dropNew)

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

  private val throttledFlow = throttle map { t =>
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

  override def onBackpressure(): Unit = {
    incrementBackpressure()
    serviceRegistryActor ! LoadMetric("IO", HighLoad)
    // Because this method will be called every time we backpressure, the timer will be overridden every
    // time until we're not backpressuring anymore
    timers.startSingleTimer(BackPressureTimerResetKey, BackPressureTimerResetAction, LoadConfig.IoNormalWindow)
  }
  
  override def actorReceive: Receive = {
    /* GCS Batch command with context */
    case (clientContext: Any, gcsBatchCommand: GcsBatchIoCommand[_, _]) =>
      val replyTo = sender()
      val commandContext= GcsBatchCommandContext(gcsBatchCommand, replyTo, Option(clientContext))
      sendToStream(commandContext)

    /* GCS Batch command without context */
    case gcsBatchCommand: GcsBatchIoCommand[_, _] =>
      val replyTo = sender()
      val commandContext= GcsBatchCommandContext(gcsBatchCommand, replyTo)
      sendToStream(commandContext)

    /* Default command with context */
    case (clientContext: Any, command: IoCommand[_]) =>
      val replyTo = sender()
      val commandContext= DefaultCommandContext(command, replyTo, Option(clientContext))
      sendToStream(commandContext)
      
    /* Default command without context */
    case command: IoCommand[_] => 
      val replyTo = sender()
      val commandContext= DefaultCommandContext(command, replyTo)
      sendToStream(commandContext)
    case BackPressureTimerResetAction => serviceRegistryActor ! LoadMetric("IO", NormalLoad)
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
  
  private val ioConfig = ConfigFactory.load().getConfig("system.io")
  
  /** Maximum number of times a command will be attempted: First attempt + 5 retries */
  val MaxAttemptsNumber: Int = ioConfig.getOrElse("number-of-attempts", 5)

  case class DefaultCommandContext[T](request: IoCommand[T], replyTo: ActorRef, override val clientContext: Option[Any] = None) extends IoCommandContext[T]
  
  case object BackPressureTimerResetKey
  case object BackPressureTimerResetAction extends ControlMessage

  def props(queueSize: Int,
            nioParallelism: Int,
            gcsParallelism: Int,
            throttle: Option[Throttle],
            serviceRegistryActor: ActorRef,
            applicationName: String,
           )
           (implicit materializer: ActorMaterializer): Props = {
    Props(new IoActor(queueSize, nioParallelism, gcsParallelism, throttle, serviceRegistryActor, applicationName)).withDispatcher(IoDispatcher)
  }
}
