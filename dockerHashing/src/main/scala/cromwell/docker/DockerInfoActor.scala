package cromwell.docker

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.stream._
import cats.effect.IO
import cats.instances.list._
import cats.syntax.parallel._
import com.google.common.cache.CacheBuilder
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.core.actor.StreamIntegration.{BackPressure, StreamContext}
import cromwell.core.{Dispatcher, DockerConfiguration}
import cromwell.docker.DockerInfoActor._
import cromwell.docker.registryv2.flows.alibabacloudcrregistry._
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import cromwell.docker.registryv2.flows.dockerhub.DockerHubRegistry
import cromwell.docker.registryv2.flows.gcr.GcrRegistry
import cromwell.docker.registryv2.flows.quay.QuayRegistry
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import fs2.Pipe
import fs2.concurrent.{NoneTerminatedQueue, Queue}
import net.ceedubs.ficus.Ficus._
import org.http4s.client.blaze.BlazeClientBuilder
import org.http4s.client.middleware.{Retry, RetryPolicy}
import org.slf4j.LoggerFactory

import scala.concurrent.duration._

final class DockerInfoActor(
                             dockerRegistryFlows: Seq[DockerRegistry],
                             queueBufferSize: Int,
                             cacheEntryTTL: FiniteDuration,
                             cacheSize: Long
                           ) extends Actor with ActorLogging {

  implicit val system = context.system
  implicit val ec = context.dispatcher
  implicit val cs = IO.contextShift(ec)
  val retryPolicy = RetryPolicy[IO](RetryPolicy.exponentialBackoff(DockerConfiguration.instance.maxTimeBetweenRetries, DockerConfiguration.instance.maxRetries))

  /* Use the guava CacheBuilder class that implements a thread safe map with built in cache features.
   * https://google.github.io/guava/releases/20.0/api/docs/com/google/common/cache/CacheBuilder.html
   * Elements will be:
   *  - Added to the cache by the streamSink, running on a thread from the stream thread pool
   *  - Accessed by this actor on its receive method thread, to check if an element is in the cache and use it
   *  - Automatically removed from the cache by the cache itself
   *  
   *  + Concurrency level is the number of expected threads to modify the cache.
   *      Set to "2" because the stream will add elements, and the cache itself remove them.
   *      This value has not a critical impact: 
   *        https://google.github.io/guava/releases/20.0/api/docs/com/google/common/cache/CacheBuilder.html#concurrencyLevel-int-
   *        
   *  + expireAfterWrite sets the time after which cache entries must expire. 
   *      We use expireAfterWrite (as opposed to expireAfterAccess because we want to the entry to expire
   *      even if it's accessed. The goal here is to force the actor to ask again for the hash after a certain 
   *      amount of time to guarantee its relative accuracy.
   *      
   *  + maximumSize sets the maximum amount of entries the cache can contain. 
   *    If/when this size is reached, least used entries will be expired
   */
  private val cache = CacheBuilder.newBuilder()
    .concurrencyLevel(2)
    .expireAfterWrite(cacheEntryTTL._1, cacheEntryTTL._2)
    .maximumSize(cacheSize)
    .build[DockerImageIdentifier, DockerInformation]()

  private def checkCache(dockerHashRequest: DockerInfoRequest) = {
    Option(cache.getIfPresent(dockerHashRequest.dockerImageID))
  }

  override def receive = receiveBehavior(Map.empty)

  def receiveBehavior(registries: Map[DockerRegistry, StreamQueue]): Receive = {
    case request: DockerInfoRequest =>
      val replyTo = sender()

      checkCache(request) match {
        case Some(cacheHit) => replyTo ! DockerInfoSuccessResponse(cacheHit, request)
        case None => sendToStream(registries, DockerInfoContext(request, replyTo))
      }
    case ShutdownCommand =>
      // Shutdown all streams by sending None to the queue
      registries
        .values.toList
        .parTraverse[IO, IO.Par, Unit](_.enqueue1(None))
        .unsafeRunSync()

      context stop self
  }

  def sendToStream(registries: Map[DockerRegistry, StreamQueue], context: DockerInfoActor.DockerInfoContext) = {
    registries collectFirst {
      case (registry, queue) if registry.accepts(context.dockerImageID) => queue
    } match {
      case Some(queue) => enqueue(context, queue)
      case None => context.replyTo ! DockerHashUnknownRegistry(context.request)
    }
  }

  def enqueue(dockerInfoContext: DockerInfoContext, queue: StreamQueue) = {
    val enqueueIO = queue.offer1(Option(dockerInfoContext)).runAsync {
      case Right(true) => IO.unit// Good !
      case _ => backpressure(dockerInfoContext)
    }

    enqueueIO.unsafeRunSync()
  }

  private def backpressure(commandContext: DockerInfoContext) = IO {
    commandContext.replyTo ! BackPressure(commandContext.request)
  }

  /*
   * Sends back responses and adds to cache in case of success. Used as the sink for each registry stream.
   */
  val streamSink: Pipe[IO, (DockerInfoResponse, DockerInfoContext), Unit] = _ evalMap {
    case (response: DockerInfoSuccessResponse, dockerInfoContext) =>
      dockerInfoContext.replyTo ! response
      IO.pure(cache.put(dockerInfoContext.dockerImageID, response.dockerInformation))
    case (response, dockerInfoContext) =>
      dockerInfoContext.replyTo ! response
      IO.pure(())
  }

  /*
   * Create a stream, a queue and start the stream.
   * The registry and queue are returned.
   * To send requests to the stream, simply enqueue a request message.
   */
  private def startAndRegisterStream(registry: DockerRegistry): IO[(DockerRegistry, StreamQueue)] = {
    implicit val timer = IO.timer(registry.config.executionContext)
    implicit val cs = IO.contextShift(registry.config.executionContext)

    /*
     * An http4s client passed to the run method of registries so they can perform the necessary requests
     * Each registry get its own client with a pool of connection.
     * Wraps the client in a retry based on the defined retry policy
     * Takes the form of a stream so it can be integrated with the actual stream of requests.
     */
    val clientStream = BlazeClientBuilder[IO](registry.config.executionContext).stream.map(Retry[IO](retryPolicy))

    Queue.boundedNoneTerminated[IO, DockerInfoContext](queueBufferSize) map { queue =>
      val source = queue.dequeue
      // If the registry imposes throttling, debounce the stream to ensure the throttling rate is respected
      val throttledSource = registry.config.throttle.map(_.delay).map(source.debounce[IO]).getOrElse(source)

      val stream = clientStream.flatMap({ client =>
        throttledSource
          // Run requests in parallel - allow nbThreads max concurrent requests - order doesn't matter
          .parEvalMapUnordered(registry.config.nbThreads)({ request => registry.run(request)(client) })
          // Send to the sink for finalization of the result
          .through(streamSink)
      })

      // Start the stream now asynchronously. It will keep running until we terminate the queue by sending None to it
      stream.compile.drain.unsafeRunAsyncAndForget()

      registry -> queue
    }
  }

  override def preStart() = {
    // Force initialization of the header constants to make sure they're valid
    locally(DockerRegistryV2Abstract)
    
    val registries =
      dockerRegistryFlows.toList
        .parTraverse(startAndRegisterStream)
        .map(_.toMap)
        // Here we need to block and wait for the registries to be ready
        .unsafeRunSync()

    context.become(receiveBehavior(registries))
    super.preStart()
  }
}

object DockerInfoActor {
  private type StreamQueue = NoneTerminatedQueue[IO, DockerInfoContext]

  val logger = LoggerFactory.getLogger("DockerRegistry")

  /* Response Messages */
  sealed trait DockerInfoResponse {
    def request: DockerInfoRequest
  }

  case class DockerSize(compressedSize: Long) {
    def toFullSize(factor: Double): Long = (compressedSize * factor).toLong
  }

  case class DockerInformation(dockerHash: DockerHashResult, dockerCompressedSize: Option[DockerSize])
  case class DockerInfoSuccessResponse(dockerInformation: DockerInformation, request: DockerInfoRequest) extends DockerInfoResponse

  sealed trait DockerHashFailureResponse extends DockerInfoResponse {
    def reason: String
  }
  case class DockerInfoFailedResponse(failure: Throwable, request: DockerInfoRequest) extends DockerHashFailureResponse {
    override val reason = s"Failed to get docker hash for ${request.dockerImageID.fullName} ${failure.getMessage}"
  }
  case class DockerHashUnknownRegistry(request: DockerInfoRequest) extends DockerHashFailureResponse {
    override val reason = s"Registry ${request.dockerImageID.host.getOrElse("<no registry>")} is not supported"
  }
  case class DockerInfoNotFound(request: DockerInfoRequest) extends DockerHashFailureResponse {
    override val reason = s"Docker image ${request.dockerImageID.fullName} not found"
  }
  case class DockerInfoUnauthorized(request: DockerInfoRequest) extends DockerHashFailureResponse {
    override val reason = s"Unauthorized to get docker hash ${request.dockerImageID.fullName}"
  }

  /* Internal ADTs */
  case class DockerInfoContext(request: DockerInfoRequest, replyTo: ActorRef) extends StreamContext {
    val dockerImageID = request.dockerImageID
    val credentials = request.credentials
  }

  private case class EnqueueResponse(result: QueueOfferResult, dockerInfoContext: DockerInfoContext)
  private case class FailedToEnqueue(failure: Throwable, dockerInfoContext: DockerInfoContext)

  def props(dockerRegistryFlows: Seq[DockerRegistry],
            queueBufferSize: Int = 100,
            cacheEntryTTL: FiniteDuration,
            cacheSize: Long) = {
    Props(new DockerInfoActor(dockerRegistryFlows, queueBufferSize, cacheEntryTTL, cacheSize)).withDispatcher(Dispatcher.IoDispatcher)
  }

  def remoteRegistriesFromConfig(config: Config): List[DockerRegistry] = {
    import cats.instances.list._
    import cats.syntax.traverse._

    val gcrConstructor = { c: DockerRegistryConfig =>
      c.copy(throttle = c.throttle.orElse(DockerConfiguration.instance.deprecatedGcrApiQueriesPer100Seconds))
      new GcrRegistry(c)
    }

    // To add a new registry, simply add it to that list
    List(
      ("dockerhub", { c: DockerRegistryConfig => new DockerHubRegistry(c) }),
      ("gcr", gcrConstructor),
      ("quay", { c: DockerRegistryConfig => new QuayRegistry(c) }),
      ("alibabacloudcr", {c: DockerRegistryConfig => new AlibabaCloudCRRegistry(c)})
    ).traverse[ErrorOr, DockerRegistry]({
      case (configPath, constructor) => DockerRegistryConfig.fromConfig(config.as[Config](configPath)).map(constructor)
    }).unsafe("Docker registry configuration")
  }
}
