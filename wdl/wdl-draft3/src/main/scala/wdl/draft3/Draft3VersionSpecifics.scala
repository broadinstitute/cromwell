package wdl.draft3

import wdl.versioning.WdlVersionSpecifics
import wom.values.{WomFloat, WomValue}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

case object Draft3VersionSpecifics extends WdlVersionSpecifics {
  override def sizeFunctionOverride(params: Seq[Try[WomValue]], sizeFunc: WomValue => Try[Double]): Try[WomFloat] = ???
}
