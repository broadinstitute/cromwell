package wes2cromwell

import java.net.URL

import scala.concurrent.Await
import scala.concurrent.duration.Duration
import akka.actor.ActorSystem
import akka.event.Logging
import akka.http.scaladsl.Http
import akka.http.scaladsl.server.Route
import com.typesafe.config.ConfigFactory
import net.ceedubs.ficus.Ficus._

object WesServer extends App with WesRunRoutes {
  val config = ConfigFactory.load()

  val port = config.as[Int]("wes2cromwell.port")
  val interface = config.as[String]("wes2cromwell.interface")

  override implicit val system: ActorSystem = ActorSystem("wes2cromwell")

  override val log = Logging(system, getClass)

  lazy val cromwellScheme = config.as[String]("cromwell.scheme")
  lazy val cromwellInterface = config.as[String]("cromwell.interface")
  lazy val cromwellPort = config.as[Int]("cromwell.port")

  override lazy val cromwellUrl = new URL(s"$cromwellScheme://$cromwellInterface:$cromwellPort")
  override val cromwellApiVersion = "v1"

  val routes: Route = runRoutes

  Http().bindAndHandle(routes, interface, port)

  println(s"Server online. Listening at port:$port")

  Await.result(system.whenTerminated, Duration.Inf)
}
