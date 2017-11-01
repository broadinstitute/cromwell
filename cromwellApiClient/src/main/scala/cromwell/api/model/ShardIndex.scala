package cromwell.api.model

import spray.json.{DefaultJsonProtocol, JsNumber, JsValue, RootJsonFormat}

case class ShardIndex(index: Option[Int]) extends AnyVal {
  override def toString: String = index.getOrElse(-1).toString
}

object ShardIndexFormatter extends DefaultJsonProtocol {
  implicit object ShardIndexJsonFormat extends RootJsonFormat[ShardIndex] {
    def write(si: ShardIndex) = JsNumber(si.index.getOrElse(-1))
    def read(value: JsValue) = value match {
      case JsNumber(i) if i.equals(-1) => ShardIndex(None)
      case JsNumber(i) if i.isValidInt && i.intValue > 0 => ShardIndex(Option(i.intValue()))
      case other => throw new UnsupportedOperationException(s"Cannot deserialize $other into a ShardIndex")
    }
  }
}
