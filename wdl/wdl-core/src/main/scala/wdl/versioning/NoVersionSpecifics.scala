package wdl.versioning

import wom.values.{WomFloat, WomValue}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

case object NoVersionSpecifics extends WdlVersionSpecifics {
  override def sizeFunctionOverride(params: Seq[Try[WomValue]], sizeFunc: WomValue => Try[Double]): Try[WomFloat] = ???
}
