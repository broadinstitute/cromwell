package cromwell.backend.io

import wom.expression.IoFunctionSet

import scala.concurrent.Future

trait FileEvaluationIoFunctionSet { this: IoFunctionSet =>
  override def glob(pattern: String) = Future.failed(new IllegalStateException("Cannot perform globing while evaluating files"))
}
