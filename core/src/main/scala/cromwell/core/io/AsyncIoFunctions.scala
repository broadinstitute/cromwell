package cromwell.core.io

import wom.expression.IoFunctionSet

import scala.concurrent.ExecutionContext

trait AsyncIoFunctions { this: IoFunctionSet =>
  /**
    * Used to perform io functions asynchronously through the ioActorProxy
    */
  def asyncIo: AsyncIo
}
