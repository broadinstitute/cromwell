package cromwell.engine.io

import java.net.{SocketException, SocketTimeoutException}
import javax.net.ssl.SSLException

import akka.NotUsed
import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.stream._
import akka.stream.scaladsl.{Flow, GraphDSL, Merge, Partition, Source}
import com.google.cloud.storage.StorageException
import com.typesafe.config.ConfigFactory
import cromwell.core.Dispatcher
import cromwell.core.actor.StreamActorHelper
import cromwell.core.actor.StreamIntegration.StreamContext
import cromwell.core.io.{IoAck, IoCommand, Throttle}
import cromwell.engine.io.IoActor._
import cromwell.engine.io.gcs.GcsBatchFlow.BatchFailedException
import cromwell.engine.io.gcs.{GcsBatchCommandContext, ParallelGcsBatchFlow}
import cromwell.engine.io.nio.NioFlow
import cromwell.filesystems.gcs.batch.GcsBatchIoCommand

/**
  * Actor that performs IO operations asynchronously using akka streams
  * 
  * @param queueSize size of the queue
  * @param throttle optional throttler to control the throughput of requests.
  *                 Applied to ALL incoming requests
  * @param materializer actor materializer to run the stream
  */
final class IoActor(queueSize: Int, throttle: Option[Throttle])(implicit val materializer: ActorMaterializer) extends Actor with ActorLogging with StreamActorHelper[IoCommandContext[_]] {
  
  implicit private val system = context.system
  
  private [io] lazy val defaultFlow = new NioFlow(parallelism = 100, context.system.scheduler).flow
  private [io] lazy val gcsBatchFlow = new ParallelGcsBatchFlow(parallelism = 10, batchSize = 100, context.system.scheduler).flow
  
  protected val source = Source.queue[IoCommandContext[_]](queueSize, OverflowStrategy.dropNew)

  protected val flow = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._
    
    val input = builder.add(Flow[IoCommandContext[_]])
    
    // Partitions requests between gcs batch, and single nio requests
    val batchPartitioner = builder.add(Partition[IoCommandContext[_]](2, {
      case gcsBatch: GcsBatchCommandContext[_, _] => 0
      case other => 1
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

  protected val throttledFlow = throttle map { t => 
    Flow[IoCommandContext[_]]
      .throttle(t.elements, t.per, t.maximumBurst, ThrottleMode.Shaping)
      .via(flow)
  } getOrElse flow
  
  override protected lazy val streamSource = source.via(throttledFlow).withAttributes(ActorAttributes.dispatcher(Dispatcher.IoDispatcher))
  
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
  }
}

trait IoCommandContext[T] extends StreamContext {
  def request: IoCommand[T]
  def replyTo: ActorRef
  def fail(failure: Throwable): IoResult = (request.fail(failure), this)
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
  val MaxAttemptsNumber = ioConfig.getOrElse[Int]("number-of-attempts", 5)

  case class DefaultCommandContext[T](request: IoCommand[T], replyTo: ActorRef, override val clientContext: Option[Any] = None) extends IoCommandContext[T]

  /**
    * ATTENTION: Transient failures are retried *forever* 
    * Be careful when adding error codes to this method.
    * Currently only 429 (= quota exceeded are considered truly transient)
    */
  def isTransient(failure: Throwable): Boolean = failure match {
    case gcs: StorageException => gcs.getCode == 429
    case _ => false
  }

  val AdditionalRetryableHttpCodes = List(
    // HTTP 410: Gone
    // From Google doc (https://cloud.google.com/storage/docs/json_api/v1/status-codes):
    // "You have attempted to use a resumable upload session that is no longer available.
    // If the reported status code was not successful and you still wish to upload the file, you must start a new session."
    410,
    // Some 503 errors seem to yield "false" on the "isRetryable" method because they are not retried.
    // The CloudStorage exception mechanism is not flawless yet (https://github.com/GoogleCloudPlatform/google-cloud-java/issues/1545)
    // so that could be the cause.
    // For now explicitly lists 503 as a retryable code here to work around that.
    503
  )
  
  /**
    * Failures that are considered retryable.
    * Retrying them should increase the "retry counter"
    */
  def isRetryable(failure: Throwable): Boolean = failure match {
    case gcs: StorageException => gcs.isRetryable || AdditionalRetryableHttpCodes.contains(gcs.getCode) || isRetryable(gcs.getCause)
    case _: SSLException => true
    case _: BatchFailedException => true
    case _: SocketException => true
    case _: SocketTimeoutException => true
    case other => isTransient(other)
  }

  def isFatal(failure: Throwable) = !isRetryable(failure)
  
  def props(queueSize: Int, throttle: Option[Throttle])(implicit materializer: ActorMaterializer) = Props(new IoActor(queueSize, throttle))
}
