package cromwell.server

import akka.actor.{ActorContext, ActorLogging, Props}
import akka.http.scaladsl.Http
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import akka.stream.ActorMaterializer
import common.util.VersionUtil
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.services.instrumentation.CromwellInstrumentationActor
import cromwell.webservice.SwaggerService
import cromwell.webservice.routes.CromwellApiService
import cromwell.webservice.routes.wes.WesRouteSupport

import scala.concurrent.Future
import scala.util.{Failure, Success}

// Note that as per the language specification, this is instantiated lazily and only used when necessary (i.e. server mode)
object CromwellServer {
  def run(gracefulShutdown: Boolean, abortJobsOnTerminate: Boolean)(cromwellSystem: CromwellSystem): Future[Any] = {
    implicit val actorSystem = cromwellSystem.actorSystem
    implicit val materializer = cromwellSystem.materializer
    actorSystem.actorOf(CromwellServerActor.props(cromwellSystem, gracefulShutdown, abortJobsOnTerminate), "cromwell-service")
    actorSystem.whenTerminated
  }
}

class CromwellServerActor(cromwellSystem: CromwellSystem, gracefulShutdown: Boolean, abortJobsOnTerminate: Boolean)(override implicit val materializer: ActorMaterializer)
  extends CromwellRootActor(
    terminator = cromwellSystem,
    gracefulShutdown = gracefulShutdown,
    abortJobsOnTerminate = abortJobsOnTerminate,
    serverMode = true,
    config = cromwellSystem.config
  )
    with CromwellApiService
    with CromwellInstrumentationActor
    with WesRouteSupport
    with SwaggerService
    with ActorLogging {
  implicit val actorSystem = context.system
  override implicit val ec = context.dispatcher
  override def actorRefFactory: ActorContext = context

  val webserviceConf = cromwellSystem.config.getConfig("webservice")
  val interface = webserviceConf.getString("interface")
  val port = webserviceConf.getInt("port")

  /**
    * /api routes have special meaning to devops' proxy servers. NOTE: the oauth mentioned on the /api endpoints in
    * cromwell.yaml is broken unless the swagger index.html is patched. Copy/paste the code from rawls or cromiam if
    * actual cromwell+swagger+oauth+/api support is needed.
    */
  val apiRoutes: Route = pathPrefix("api")(concat(workflowRoutes, womtoolRoutes))
  val nonApiRoutes: Route = concat(engineRoutes, swaggerUiResourceRoute, wesRoutes)
  val allRoutes: Route = concat(apiRoutes, nonApiRoutes)

  val serverBinding = Http().bindAndHandle(allRoutes, interface, port)

  CromwellShutdown.registerUnbindTask(actorSystem, serverBinding)

  serverBinding onComplete {
    case Success(serverBindingActual) =>
      val version = VersionUtil.getVersion("cromwell-engine")
      val host = serverBindingActual.localAddress.getHostString
      val port = serverBindingActual.localAddress.getPort
      actorSystem.log.info(s"Cromwell {} service started on {}:{}...", version, host, port)
    case Failure(e) =>
      /*
        TODO:
        If/when CromwellServer behaves like a better async citizen, we may be less paranoid about our async log messages
        not appearing due to the actor system shutdown. For now, synchronously print to the stderr so that the user has
        some idea of why the server failed to start up.
      */
      Console.err.println(s"Binding failed interface $interface port $port")
      e.printStackTrace(Console.err)
      cromwellSystem.shutdownActorSystem()
  }

  /*
    During testing it looked like not explicitly invoking the WMA in order to evaluate all of the lazy actors in
    CromwellRootActor would lead to weird timing issues the first time it was invoked organically
   */
  workflowManagerActor
}

object CromwellServerActor {
  def props(cromwellSystem: CromwellSystem, gracefulShutdown: Boolean, abortJobsOnTerminate: Boolean)(implicit materializer: ActorMaterializer): Props = {
    Props(new CromwellServerActor(cromwellSystem, gracefulShutdown, abortJobsOnTerminate)).withDispatcher(EngineDispatcher)
  }
}
