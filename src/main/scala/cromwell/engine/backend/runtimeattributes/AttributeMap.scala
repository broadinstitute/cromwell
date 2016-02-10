package cromwell.engine.backend.runtimeattributes

import AttributeMap.EnhancedBackendType
import cromwell.engine.backend.BackendType

object AttributeMap {
  implicit class EnhancedBackendType(val backendType: BackendType) extends AnyVal {
    def supportedKeys: Set[RuntimeKey] = for {
      key <- RuntimeKey.values().toSet
      if key.supports(backendType)
    } yield key
  }
}

trait AttributeMapTrait {
  def keys: Set[String]
  def get(key: RuntimeKey): Option[String]
  def getSeq(key: RuntimeKey): Option[Seq[String]]
  def unsupportedKeys(backendType: BackendType): Seq[String]
}

case class AttributeMap(attrs: Map[String, Seq[String]]) extends AttributeMapTrait {
  def keys = attrs.keySet
  def get(key: RuntimeKey): Option[String] = getSeq(key).flatMap(_.headOption)
  def getSeq(key: RuntimeKey): Option[Seq[String]] = attrs.get(key.key)

  def unsupportedKeys(backendType: BackendType): Seq[String] = {
    val supportedKeys = backendType.supportedKeys map { _.key }
    val unsupportedKeys = attrs.keySet -- supportedKeys

    if (unsupportedKeys.isEmpty) Vector.empty
    else Vector(s"Found unsupported keys for backend '$backendType': " + unsupportedKeys.toSeq.sorted.mkString(", "))
  }
}


