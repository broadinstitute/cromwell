package cromwell.docker.registryv2.flows

import akka.stream.FanOutShape2
import akka.stream.scaladsl.{GraphDSL, Partition}

import scala.util.{Failure, Success, Try}

object FlowUtils {

  /**
    * Takes an input of the form (Try[T], U) and exposes 2 output ports.
    * U is the type of the context value to be passed along.
    * The first one will emit a pair of the form (value, context) if the try is a success.
    * The second one will emit a pair of the form (throwable, context) if the try is a failure.
    */
  def fanOutTry[T, U] = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._

    val partition = builder.add(Partition[(Try[T], U)](2, {
      case (Success(_), _) => 0
      case (Failure(_), _) => 1
    }))

    val successOut: PortOps[(T, U)] = partition.out(0) collect {
      case (Success(value), flowContext) => (value, flowContext)
    }

    val failureOut: PortOps[(Throwable, U)] = partition.out(1) collect {
      case (Failure(failure), flowContext) => (failure, flowContext)
    }

    new FanOutShape2[(Try[T], U), (T, U), (Throwable, U)](partition.in, successOut.outlet, failureOut.outlet)
  }
}
