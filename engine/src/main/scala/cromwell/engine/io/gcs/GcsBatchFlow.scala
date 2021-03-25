package cromwell.engine.io.gcs

import akka.NotUsed

import java.io.IOException
import akka.actor.Scheduler
import akka.stream._
import akka.stream.scaladsl.{Flow, GraphDSL, MergePreferred, Partition}
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.http.{HttpRequest, HttpRequestInitializer}
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.storage.Storage
import com.typesafe.scalalogging.StrictLogging
import common.util.StringUtil.EnhancedToStringable
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.core.io.IoAck
import cromwell.engine.io.IoActor._
import cromwell.engine.io.IoAttempts.EnhancedCromwellIoException
import cromwell.engine.io.RetryableRequestSupport.{isInfinitelyRetryable, isRetryable}
import cromwell.engine.io.gcs.GcsBatchFlow._
import cromwell.engine.io.{IoAttempts, IoCommandContext}
import mouse.boolean._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Try}

object GcsBatchFlow {

  /**
    * Exception used to fail the request promises when the batch request itself fails.
    * Is considered retryable.
    */
  case class BatchFailedException(failure: Throwable) extends IOException(failure)

  private val ReadForbiddenPattern = ".*does not have storage\\.objects\\.(?:get|list|copy) access to ([^/]+).*".r.pattern

  /* Returns `Some(bucket)` if the specified argument represents a forbidden attempt to read from `bucket`. */
  private[gcs] def getReadForbiddenBucket(errorMsg: String): Option[String] = {
    val matcher = ReadForbiddenPattern.matcher(errorMsg)
    matcher.matches().option(matcher.group(1))
  }
}

class GcsBatchFlow(batchSize: Int, scheduler: Scheduler, onRetry: IoCommandContext[_] => Throwable => Unit, applicationName: String)
                  (implicit ec: ExecutionContext) extends StrictLogging {

  // Does not carry any authentication, assumes all underlying requests are properly authenticated
  private val httpRequestInitializer = new HttpRequestInitializer {
    override def initialize(request: HttpRequest): Unit = {
      request.setConnectTimeout(GoogleConfiguration.DefaultConnectionTimeout.toMillis.toInt)
      request.setReadTimeout(GoogleConfiguration.DefaultReadTimeout.toMillis.toInt)
      ()
    }
  }

  /**
    * Returns a new BatchRequest instance.
    *
    * The correlated errors that were seen were:
    * - Via CROM-6708 timeouts were occurring, followed usually by...
    * - Via CROM-6709 batches had grown beyond 3001+ requests.
    *
    * Examining BatchRequest.execute(), if any exception occurs during the method, either during the HTTP request or
    * during one of the response handlers, the BatchRequest.execute() method simply throws the exception.
    *
    * When exceptions are NOT thrown the BatchRequest's internal queue is only partially or fully cleared:
    * - partial : https://github.com/googleapis/google-api-java-client/blob/v1.31.1/google-api-client/src/main/java/com/google/api/client/googleapis/batch/BatchRequest.java#L275
    * -   full  : https://github.com/googleapis/google-api-java-client/blob/v1.31.1/google-api-client/src/main/java/com/google/api/client/googleapis/batch/BatchRequest.java#L281
    *
    * The BatchRequest's internal queue is NOT cleared when exceptions are thrown inside BatchRequest.execute().
    * Any subsequent enqueueing operations only add to the internal queue and do not replace the internal queue.
    * So eventually timeouts lead to large batches of 1000+ elements that cause further timeout exceptions too.
    * Those failed batch elements were re-appended leading to batches of greater than 3000 requests.
    * At least that's the theory.
    *
    * Instead we'll create a new BatchRequest each time we run GcsBatchFlow.executeBatch().
    *
    * From basic performance testing with YourKit we were able to create 1,000,000 of these BatchRequest objects in
    * under 17s. However, if needed a cached var batchRequest instance could be implemented that only recreates the
    * object when the internal queue is "dirty" with a .size() > 0.
    */
  private def newBatchRequest(): BatchRequest = {
    val builder = new Storage.Builder(
      GcsStorage.HttpTransport,
      JacksonFactory.getDefaultInstance,
      httpRequestInitializer
    ).setApplicationName(applicationName)
    val client = builder.build()
    client.batch(client.getRequestFactory.getInitializer)
  }

  val flow: Graph[
    FlowShape[
      GcsBatchCommandContext[_, _],
      (IoAck[_], IoCommandContext[_]),
    ],
    NotUsed,
  ] = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._

    // Source where batch commands are coming from. This is the input port of this flow
    val source = builder.add(Flow[GcsBatchCommandContext[_, _]])

    // Merge commands from source (above), and commands that need to be retried (see retries below)
    val sourceMerger = builder.add(MergePreferred[GcsBatchCommandContext[_, _]](1))

    // Process a batch and spit atomic GcsBatchResponses out for each internal request
    val batchProcessor = builder.add(
      Flow[GcsBatchCommandContext[_, _]]
        // Group commands together in batches so they can be processed as such
      .groupedWithin(batchSize, 5 seconds)
        // execute the batch and outputs each sub-response individually, as a Future
      .mapConcat[Future[GcsBatchResponse[_]]](executeBatch)
        // Wait for each Future to complete
      .mapAsyncUnordered[GcsBatchResponse[_]](batchSize) { identity }
    )

    // Partitions the responses: Terminal responses exit the flow, others go back to the sourceMerger
    val responseHandler = builder.add(responseHandlerFlow)

    // Buffer commands to be retried to avoid backpressuring too rapidly
    val nextRequestBuffer = builder.add(Flow[GcsBatchCommandContext[_, _]].buffer(batchSize, OverflowStrategy.backpressure))

    source ~> sourceMerger ~> batchProcessor ~> responseHandler.in
              sourceMerger.preferred <~ nextRequestBuffer <~ responseHandler.out1

    FlowShape[GcsBatchCommandContext[_, _], IoResult](source.in, responseHandler.out0)
  }

  /**
    * Fan out shape splitting GcsBatchResponse into 2:
    *   First port emits terminal result that can exit the GcsBatch flow
    *   Second port emits request to be re-injected to be executed in a later batch
    */
  private lazy val responseHandlerFlow = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._

    val source = builder.add(Partition[GcsBatchResponse[_]](2, {
      case _: GcsBatchTerminal[_] => 0
      case _ => 1
    }))

    // Terminal responses: output of this flow
    val terminals = source.out(0) collect { case terminal: GcsBatchTerminal[_] => terminal.ioResult }

    // Next command context, can be a retry or another request needed by the command
    val nextRequest = source.out(1).collect {
      case retry: GcsBatchRetry[_] => retry.context
      case nextRequest: GcsBatchNextRequest[_] => nextRequest.context
    }

    new FanOutShape2[GcsBatchResponse[_], IoResult, GcsBatchCommandContext[_, _]](source.in, terminals.outlet, nextRequest.outlet)
  }

  private def executeBatch(contexts: Seq[GcsBatchCommandContext[_, _]]): List[Future[GcsBatchResponse[_]]] = {
    def failAllPromisesWith(failure: Throwable): Unit = contexts foreach { context =>
      context.promise.tryFailure(failure)
      ()
    }

    val batchRequest = newBatchRequest()

    // Add all requests to the batch
    contexts foreach { _.queue(batchRequest) }

    val batchCommandNamesList = contexts.map(_.request.toString)
    // Try to execute the batch request.
    // If it fails with an IO Exception, fail all the underlying promises with a retryable BatchFailedException
    // Otherwise fail with the original exception
    Try(batchRequest.execute()) match {
      case Failure(failure: IOException) =>
        logger.info(s"Failed to execute GCS Batch request. Failed request belonged to batch of size ${batchRequest.size()} containing commands: " +
          s"${batchCommandNamesList.mkString("\n")}.\n${failure.toPrettyElidedString(limit = 1000)}")
        failAllPromisesWith(BatchFailedException(failure))
      case Failure(failure) =>
        logger.info(s"Failed to execute GCS Batch request. Failed request belonged to batch of size ${batchRequest.size()} containing commands: " +
          s"${batchCommandNamesList.mkString("\n")}.\n${failure.toPrettyElidedString(limit = 1000)}")
        failAllPromisesWith(failure)
      case _ =>
    }

    // Map all promise responses to a GcsBatchResponse to be either sent back as a response or retried in the next batch
    contexts.toList map { context =>
      context.promise.future map {
        case Left(response) => GcsBatchTerminal(response)
        case Right(nextRequest) => GcsBatchNextRequest(nextRequest)
      } recoverWith recoverCommand(context)
    }
  }

  /**
    * Handles a failed future.
    *   If the failure is retryable, and the command hasn't reached its max attempts:
    *     schedule the command to be retried in a later batch after waiting for the appropriate amount of time
    *   Otherwise create a GcsBatchTerminal response with the IoFailure
    *   In both cases, returns a successful Future to avoid failing the stream or dropping elements
    */
  private def recoverCommand(context: GcsBatchCommandContext[_, _]): PartialFunction[Throwable, Future[GcsBatchResponse[_]]] = {
    // If the failure is retryable - recover with a GcsBatchRetry so it can be retried in the next batch
    case failure if isRetryable(failure) =>
      context.retryIn match {
        case Some(waitTime) if isInfinitelyRetryable(failure) =>
          onRetry(context)(failure)
          akka.pattern.after(waitTime, scheduler)(Future.successful(GcsBatchRetry(context.nextTransient, failure)))
        case Some(waitTime) =>
          onRetry(context)(failure)
          akka.pattern.after(waitTime, scheduler)(Future.successful(GcsBatchRetry(context.next, failure)))
        case None => fail(context, failure)
      }
    // Otherwise just fail the command, either with a specific "read forbidden" failure or just generic failure.
    case failure =>
      Option(failure.getMessage) match {
        case Some(errorMsg) =>
          getReadForbiddenBucket(errorMsg) match {
            case Some(bucket) => failReadForbidden(context, failure, bucket)
            case None => fail(context, failure)
          }
        case None => fail(context, failure)
      }
  }

  /**
    * Fail a command context with a failure.
    */
  private def fail(context: GcsBatchCommandContext[_, _], failure: Throwable) = {
    Future.successful(
      GcsBatchTerminal(
        context.fail(EnhancedCromwellIoException(IoAttempts(context.currentAttempt), failure))
      )
    )
  }

  /**
    * Fail a command context with a forbidden failure.
    */
  private def failReadForbidden(context: GcsBatchCommandContext[_, _], failure: Throwable, forbiddenPath: String) = {
    Future.successful(
      GcsBatchTerminal(
        context.failReadForbidden(EnhancedCromwellIoException(IoAttempts(context.currentAttempt), failure), forbiddenPath)
      )
    )
  }
}
