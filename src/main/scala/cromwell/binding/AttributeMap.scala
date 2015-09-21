package cromwell.binding

import cromwell.parser.{BackendType, RuntimeKey}
import AttributeMap._

object AttributeMap {
  implicit class EnhancedBackendType(val backendType: BackendType) extends AnyVal {
    def supportedKeys: Set[RuntimeKey] = for {
      key <- RuntimeKey.values().toSet
      if key.supports(backendType)
    } yield key
  }
}

case class AttributeMap(attrs: Map[String, Seq[String]]) {
  def get(key: RuntimeKey): Option[String] = attrs.get(key.key).flatMap(_.headOption)
  def getSeq(key: RuntimeKey): Option[Seq[String]] = attrs.get(key.key)

  def unsupportedKeys(backendType: BackendType): Seq[String] = {
    val supportedKeys = backendType.supportedKeys map { _.key }
    val unsupportedKeys = attrs.keySet -- supportedKeys

    if (unsupportedKeys.isEmpty) Vector.empty
    else Vector(s"Found unsupported keys for backend '$backendType': " + unsupportedKeys.toSeq.sorted.mkString(", "))
  }
}


