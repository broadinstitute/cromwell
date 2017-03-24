package cromwell.server

import java.util.concurrent.TimeoutException

import akka.actor.{ActorContext, ActorSystem, Props}
import com.typesafe.config.Config
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.webservice.WorkflowJsonSupport._
import cromwell.webservice.{APIResponse, CromwellApiService, SwaggerService}
import cromwell.webservice.SprayCanHttpService._
import cromwell.webservice.WrappedRoute._
import net.ceedubs.ficus.Ficus._
import spray.http._
import spray.json._
import spray.routing.Route

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContextExecutor, Future}
import scala.util.{Failure, Success}

// Note that as per the language specification, this is instantiated lazily and only used when necessary (i.e. server mode)
object CromwellServer {

  def run(cromwellSystem: CromwellSystem): Future[Any] = {
    implicit val executionContext = scala.concurrent.ExecutionContext.Implicits.global

    val actorSystem: ActorSystem = cromwellSystem.actorSystem

    val service = actorSystem.actorOf(CromwellServerActor.props(cromwellSystem.conf), "cromwell-service")
    val webserviceConf = cromwellSystem.conf.getConfig("webservice")

    val interface = webserviceConf.getString("interface")
    val port = webserviceConf.getInt("port")
    val timeout = webserviceConf.as[FiniteDuration]("binding-timeout")
    val futureBind = service.bind(interface, port)(implicitly, timeout, actorSystem, implicitly)
    futureBind andThen {
      case Success(_) =>
        actorSystem.log.info("Cromwell service started...")
        Await.result(actorSystem.whenTerminated, Duration.Inf)
      case Failure(throwable) =>
        /*
        TODO:
        If/when CromwellServer behaves like a better async citizen, we may be less paranoid about our async log messages
        not appearing due to the actor system shutdown. For now, synchronously print to the stderr so that the user has
        some idea of why the server failed to start up.
         */
        Console.err.println(s"Binding failed interface $interface port $port")
        throwable.printStackTrace(Console.err)
        cromwellSystem.shutdownActorSystem()
    }
  }
}

class CromwellServerActor(config: Config) extends CromwellRootActor with CromwellApiService with SwaggerService {
  implicit def executionContext: ExecutionContextExecutor = actorRefFactory.dispatcher

  override val serverMode = true
  override val abortJobsOnTerminate = false

  override def actorRefFactory: ActorContext = context
  override def receive: PartialFunction[Any, Unit] = handleTimeouts orElse runRoute(possibleRoutes)

  val routeUnwrapped: Boolean = config.as[Option[Boolean]]("api.routeUnwrapped").getOrElse(false)
  val possibleRoutes: Route = workflowRoutes.wrapped("api", routeUnwrapped) ~ swaggerUiResourceRoute
  val timeoutError: String = APIResponse.error(new TimeoutException(
    "The server was not able to produce a timely response to your request.")).toJson.prettyPrint

  def handleTimeouts: Receive = {
    case Timedout(_: HttpRequest) =>
      sender() ! HttpResponse(StatusCodes.InternalServerError, HttpEntity(ContentType(MediaTypes.`application/json`), timeoutError))
  }

  /*
    During testing it looked like not explicitly invoking the WMA in order to evaluate all of the lazy actors in
    CromwellRootActor would lead to weird timing issues the first time it was invoked organically
   */
  workflowManagerActor
}

object CromwellServerActor {
  def props(config: Config): Props = {
    Props(new CromwellServerActor(config)).withDispatcher(EngineDispatcher)
  }
}
