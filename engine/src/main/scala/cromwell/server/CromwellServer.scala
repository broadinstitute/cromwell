package cromwell.server

import akka.util.Timeout
import cromwell.webservice.CromwellApiServiceActor
import lenthall.spray.SprayCanHttpService._

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.util.{Failure, Success}

// Note that as per the language specification, this is instantiated lazily and only used when necessary (i.e. server mode)
object CromwellServer {
  implicit val timeout = Timeout(5.seconds)
  import scala.concurrent.ExecutionContext.Implicits.global


  def run(workflowManagerSystem: WorkflowManagerSystem): Future[Any] = {
    implicit val actorSystem = workflowManagerSystem.actorSystem

    val service = actorSystem.actorOf(CromwellApiServiceActor.props(
      workflowManagerSystem.workflowManagerActor,
      workflowManagerSystem.workflowStoreActor,
      workflowManagerSystem.conf), "cromwell-service")
    val webserviceConf = workflowManagerSystem.conf.getConfig("webservice")

    val interface = webserviceConf.getString("interface")
    val port = webserviceConf.getInt("port")
    val futureBind = service.bind(interface = interface, port = port)
    futureBind andThen {
      case Success(_) =>
        actorSystem.log.info("Cromwell service started...")
        actorSystem.awaitTermination()
      case Failure(throwable) =>
        /*
        TODO:
        If/when CromwellServer behaves like a better async citizen, we may be less paranoid about our async log messages
        not appearing due to the actor system shutdown. For now, synchronously print to the stderr so that the user has
        some idea of why the server failed to start up.
         */
        Console.err.println(s"Binding failed interface $interface port $port")
        throwable.printStackTrace(Console.err)
        workflowManagerSystem.shutdownActorSystem()
    }
  }
}
