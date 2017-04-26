package cromwell.docker

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.stream._
import akka.stream.scaladsl.{GraphDSL, Merge, Partition, Sink, Source}
import com.google.common.cache.CacheBuilder
import cromwell.core.Dispatcher
import cromwell.core.actor.StreamActorHelper
import cromwell.core.actor.StreamIntegration.StreamContext
import cromwell.docker.DockerHashActor._
import org.slf4j.LoggerFactory

import scala.concurrent.duration.FiniteDuration

final class DockerHashActor(
                           dockerRegistryFlows: Seq[DockerFlow],
                           queueBufferSize: Int,
                           cacheEntryTTL: FiniteDuration,
                           cacheSize: Long
                          )(implicit val materializer: ActorMaterializer) extends Actor with ActorLogging with StreamActorHelper[DockerHashContext] {

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
    .build[DockerImageIdentifierWithoutHash, DockerHashResult]()

  /*
   * Intermediate sink responsible for updating the cache as soon as a successful hash is retrieved
   */
  private val updateCacheSink = Sink.foreach[(DockerHashResponse, DockerHashContext)] {
    case (response: DockerHashResponseSuccess, dockerHashContext) => 
      cache.put(dockerHashContext.dockerImageID, response.dockerHash)
    case _ => // Only put successful hashes in the cache
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

    val noMatch = partitionFlows.out(dockerRegistryFlows.size).map(dockerContext => (DockerHashUnknownRegistry(dockerContext.request), dockerContext)).outlet
    
    noMatch ~> mergeFlows.in(dockerRegistryFlows.size)
    
    FlowShape(partitionFlows.in, mergeFlows.out)
  }

  private def checkCache(dockerHashRequest: DockerHashRequest) = {
    Option(cache.getIfPresent(dockerHashRequest.dockerImageID)) map { hashResult => 
      DockerHashResponseSuccess(hashResult, dockerHashRequest) 
    }
  }

  override protected def actorReceive: Receive = {
    case request: DockerHashRequest =>
      val replyTo = sender()

      checkCache(request) match {
        case Some(cacheHit) => replyTo ! cacheHit
        case None => sendToStream(DockerHashContext(request, replyTo))
      }
  }

  override protected def streamSource = Source.queue[DockerHashContext](queueBufferSize, OverflowStrategy.dropNew)
    .via(dockerFlow)
    .alsoTo(updateCacheSink)
    .withAttributes(ActorAttributes.dispatcher(Dispatcher.IoDispatcher))
}

object DockerHashActor {

  val logger = LoggerFactory.getLogger("DockerRegistry")

  /* Response Messages */
  sealed trait DockerHashResponse {
    def request: DockerHashRequest
  }
  
  case class DockerHashResponseSuccess(dockerHash: DockerHashResult, request: DockerHashRequest) extends DockerHashResponse
  
  sealed trait DockerHashFailureResponse extends DockerHashResponse {
    def reason: String
  }
  case class DockerHashFailedResponse(failure: Throwable, request: DockerHashRequest) extends DockerHashFailureResponse {
    override val reason = s"Failed to get docker hash for ${request.dockerImageID.fullName} ${failure.getMessage}"
  }
  case class DockerHashUnknownRegistry(request: DockerHashRequest) extends DockerHashFailureResponse {
    override val reason = s"Registry ${request.dockerImageID.host} is not supported"
  }
  case class DockerHashNotFound(request: DockerHashRequest) extends DockerHashFailureResponse {
    override val reason = s"Docker image ${request.dockerImageID.fullName} not found"
  }
  case class DockerHashUnauthorized(request: DockerHashRequest) extends DockerHashFailureResponse {
    override val reason = s"Unauthorized to get docker hash ${request.dockerImageID.fullName}"
  }

  /* Internal ADTs */
  case class DockerHashContext(request: DockerHashRequest, replyTo: ActorRef) extends StreamContext {
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
