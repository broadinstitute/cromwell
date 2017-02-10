package cromwell.core.callcaching.docker

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.pattern.pipe
import akka.stream.QueueOfferResult.{Dropped, Enqueued, QueueClosed}
import akka.stream._
import akka.stream.scaladsl.{GraphDSL, Merge, Partition, Sink, Source}
import com.google.common.cache.CacheBuilder
import cromwell.core.callcaching.docker.DockerHashActor._
import org.slf4j.LoggerFactory

import scala.concurrent.Future
import scala.concurrent.duration.FiniteDuration

final class DockerHashActor(
                           dockerRegistryFlows: Seq[DockerFlow],
                           queueBufferSize: Int,
                           cacheEntryTTL: FiniteDuration,
                           cacheSize: Long
                          )(implicit val materializer: ActorMaterializer) extends Actor with ActorLogging {

  implicit val system = context.system
  implicit val ec = context.dispatcher
  
  /* Use the guava CacheBuilder class that implements a thread safe map with built in cache features.
   * https://google.github.io/guava/releases/20.0/api/docs/com/google/common/cache/CacheBuilder.html
   * Elements will be:
   *  - Added to the cache by the updateCacheSink, running on a thread from the stream thread pool
   *  - Accessed by this actor on its receive method thread, to check if an element is in the cache and use it
   *  - Automatically removed from the cache by the cache itself
   *  
   *  + Concurrency level is the number of expected threads to modify the cache.
   *      Choose 2 because the stream will add elements, and the cache itself remove them.
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
    .build[DockerImageIdentifierWithoutHash, DockerHashResult]()
  
  /*
   * Never fail the stream, if an exception occurs just drop the element
   * The sender is responsible for sending a new message if it didn't get a response
   * within a reasonable amount of time
   */
  private val decider: Supervision.Decider = _ => Supervision.Resume

  /*
   * Intermediate sink responsible for updating the cache as soon as a successful hash is retrieved
   */
  private val updateCacheSink = Sink.foreach[(DockerHashResponse, DockerHashContext)] {
    case (response: DockerHashResponseSuccess, dockerHashContext) => 
      cache.put(dockerHashContext.dockerImageID, response.dockerHash)
    case _ => // Only put successful hashes in the cache
  }
  
  /*
   * Final sink of the stream. Sends the response to the caller.
   */
  private val replySink = Sink.foreach[(DockerHashResponse, DockerHashContext)] {
    case (response, dockerHashContext) => dockerHashContext.replyTo ! response
  }

  /*
    Send each message to the corresponding docker registry flow if it exists,
    and then merge the results back together.
    If no flow is able to handle the request, create a DockerHashUnknownRegistry response
                  _______                       ______
                  |      |                      |     |
                  |      | -> dockerhub flow -> |     | 
      hash        |      | -> gcr flow  ->      |     | -> Response
      context ->  |      | -> ...       ->      |     |
                  |______|                      |_____|
   */
  private val dockerFlow = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._
    
    def mapContextToFlow(dockerHashContext: DockerHashContext) = {
      dockerRegistryFlows.indexWhere(_.accepts(dockerHashContext.dockerImageID)) match {
          // requests with no matching host are sent to the last port
        case -1 => dockerRegistryFlows.length
        case found => found
      }
    }

    val partitionFlows = builder.add(Partition[DockerHashContext](dockerRegistryFlows.size + 1, mapContextToFlow))
    val mergeFlows = builder.add(Merge[(DockerHashResponse, DockerHashContext)](dockerRegistryFlows.size + 1))
    
    // Build a partial flow for each dockerRegistryFlow available
    val flows = dockerRegistryFlows map { _.buildFlow() }

    // Connect each flow to one output port of the broadcast, and merge them all back together
    dockerRegistryFlows.indices foreach { i =>
      partitionFlows.out(i) ~> flows(i) ~> mergeFlows.in(i)
    }

    val noMatch = partitionFlows.out(dockerRegistryFlows.size).map(dockerContext => (DockerHashUnknownRegistry(dockerContext.dockerImageID), dockerContext)).outlet
    
    noMatch ~> mergeFlows.in(dockerRegistryFlows.size)
    
    FlowShape(partitionFlows.in, mergeFlows.out)
  }

  private val queue = Source.queue[DockerHashContext](queueBufferSize, OverflowStrategy.dropNew)
    .via(dockerFlow)
    .alsoTo(updateCacheSink)
    .to(replySink)
    .withAttributes(ActorAttributes.supervisionStrategy(decider))
    .run()
  
  override def receive: Receive = {
    // When receiving a request, just enqueue it and pipe the result back to self
    // Note: could be worth having another actor to send it to, to prevent overflowing this mailbox
    case request: DockerHashRequest =>
      val replyTo = sender()
      
      checkCache(request.dockerImageID) match {
        case Some(cacheHit) => replyTo ! cacheHit
        case None => sendToStream(DockerHashContext(request, replyTo))
      }
    case EnqueueResponse(Enqueued, dockerHashContext) => // All good !
    case EnqueueResponse(Dropped, dockerHashContext) => dockerHashContext.replyTo ! DockerHashBackPressure(dockerHashContext.request)
      // In any of the below case, the stream is in a state where it will not be able to receive new elements.
      // This means something blew off so we can just restart the actor to re-instantiate a new stream
    case EnqueueResponse(QueueClosed, dockerHashContext) =>
      val exception = new RuntimeException(s"Failed to enqueue docker hash request ${dockerHashContext.request}. Queue was closed")
      logAndRestart(exception)
    case EnqueueResponse(QueueOfferResult.Failure(failure), dockerHashContext) => 
      logAndRestart(failure)
    case FailedToEnqueue(throwable, dockerHashContext) =>
      logAndRestart(throwable)
  }
  
  private def logAndRestart(throwable: Throwable) = {
    log.error("Failed to process docker hash request", throwable)
    // Throw the exception that will be caught by supervisor and restart the actor
    throw DockerHashActorException(throwable)
  }
  
  private def checkCache(dockerImageIdentifierWithoutHash: DockerImageIdentifierWithoutHash) = {
    Option(cache.getIfPresent(dockerImageIdentifierWithoutHash)) map { hashResult => 
      DockerHashResponseSuccess(hashResult) 
    }
  }
  
  private def sendToStream(dockerHashContext: DockerHashContext) = {
    val enqueue = queue offer dockerHashContext map { result =>
      EnqueueResponse(result, dockerHashContext)
    } recoverWith {
      case t => Future.successful(FailedToEnqueue(t, dockerHashContext))
    }

    pipe(enqueue) to self
    ()
  }
}

object DockerHashActor {

  val logger = LoggerFactory.getLogger("DockerRegistry")

  /* Response Messages */
  sealed trait DockerHashResponse
  
  case class DockerHashResponseSuccess(dockerHash: DockerHashResult) extends DockerHashResponse
  
  sealed trait DockerHashFailureResponse extends DockerHashResponse {
    def reason: String
  }
  case class DockerHashFailedResponse(failure: Throwable, dockerIdentifier: DockerImageIdentifierWithoutHash) extends DockerHashFailureResponse {
    override val reason = s"Failed to get docker hash for ${dockerIdentifier.fullName} ${failure.getMessage}"
  }
  case class DockerHashUnknownRegistry(dockerImageId: DockerImageIdentifierWithoutHash) extends DockerHashFailureResponse {
    override val reason = s"Registry ${dockerImageId.host} is not supported"
  }
  case class DockerHashNotFound(dockerImageId: DockerImageIdentifierWithoutHash) extends DockerHashFailureResponse {
    override val reason = s"Docker image ${dockerImageId.fullName} not found"
  }
  case class DockerHashUnauthorized(dockerImageId: DockerImageIdentifierWithoutHash) extends DockerHashFailureResponse {
    override val reason = s"Unauthorized to get docker hash ${dockerImageId.fullName}"
  }

  case class DockerHashBackPressure(originalRequest: DockerHashRequest) extends DockerHashResponse
  case class DockerHashActorException(failure: Throwable) extends RuntimeException(failure)
  
  /* Internal ADTs */
  case class DockerHashContext(request: DockerHashRequest, replyTo: ActorRef) {
    val dockerImageID = request.dockerImageID
    val credentials = request.credentials
  }

  private case class EnqueueResponse(result: QueueOfferResult, dockerHashContext: DockerHashContext)
  private case class FailedToEnqueue(failure: Throwable, dockerHashContext: DockerHashContext)

  def props(dockerRegistryFlows: Seq[DockerFlow],
            queueBufferSize: Int = 100,
            cacheEntryTTL: FiniteDuration,
            cacheSize: Long)(materializer: ActorMaterializer) = Props(new DockerHashActor(dockerRegistryFlows, queueBufferSize, cacheEntryTTL, cacheSize)(materializer))
}
