package cromwell.engine

import akka.util.Timeout

import scala.concurrent.duration._
import scala.language.postfixOps

trait CromwellActor {
  protected implicit val timeout = Timeout(5 seconds)
}
