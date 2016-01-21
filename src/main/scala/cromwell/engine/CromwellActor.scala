package cromwell.engine

import akka.util.Timeout
import com.typesafe.config.{ConfigException, ConfigFactory}

import scala.concurrent.duration._
import scala.language.postfixOps

trait CromwellActor {
  protected implicit val timeout = Timeout(5 seconds)

  /**
   * Retrieves the configuration option that determines whether this actor should abort all jobs if it receives
   * a shutdown hook.
   * @return - The value of the configuration option, or 'false' if the option isn't specified.
   */
  def getAbortJobsOnTerminate: Boolean = {
    val config=ConfigFactory.load.getConfig("backend")
    try {
      config.getBoolean("abortJobsOnTerminate")
    } catch {
      case _:ConfigException  => false
    }
  }

}
