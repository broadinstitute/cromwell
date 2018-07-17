package wes2cromwell

import scala.concurrent.Await
import scala.concurrent.duration.Duration
import akka.actor.{ActorRef, ActorSystem}
import akka.event.Logging
import akka.http.scaladsl.Http
import akka.http.scaladsl.server.Route
import akka.stream.ActorMaterializer
import com.typesafe.config.ConfigFactory
import net.ceedubs.ficus.Ficus._

// MAIN
object WesServer extends App with WesWorkflowRoutes {
  val config = ConfigFactory.load()

  val port = config.as[Int]("wes2cromwell.port")
  val interface = config.as[String]("wes2cromwell.interface")

  // set up ActorSystem and other dependencies here
  implicit val system: ActorSystem = ActorSystem("helloAkkaHttpServer")
  implicit val materializer: ActorMaterializer = ActorMaterializer()
  val workflowActor: ActorRef = system.actorOf(WorkflowActor.props, "workflowRegistryActor")

  override val log = Logging(system, getClass)

  lazy val cromwellScheme = config.as[String]("cromwell.scheme")
  lazy val cromwellInterface = config.as[String]("cromwell.interface")
  lazy val cromwellPort = config.as[Int]("cromwell.port")
  override lazy val cromwellPath = s"$cromwellScheme://$cromwellInterface:$cromwellPort/api/workflows/v1"

  // from the UserRoutes trait
  val routes: Route = workflowRoutes

  Http().bindAndHandle(routes, interface, port)

  println(s"Server online. Listening at port:$port")

  Await.result(system.whenTerminated, Duration.Inf)
}
