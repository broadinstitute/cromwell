package wdl.versioning

import wom.values.{WomFloat, WomValue}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

trait WdlVersionSpecifics {
  def sizeFunctionOverride(params: Seq[Try[WomValue]], sizeFunc: WomValue => Try[Double]): Try[WomFloat]
}
