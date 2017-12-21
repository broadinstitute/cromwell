package cromwell.core.io

import wom.expression.IoFunctionSet

import scala.concurrent.ExecutionContext

trait AsyncIoFunctions { this: IoFunctionSet =>
  /**
    * Used to perform io functions asynchronously through the IoActorEndpoint
    */
  def asyncIo: AsyncIo

  /**
    * To map/flatMap over IO results
    */
  implicit def ec: ExecutionContext
}
