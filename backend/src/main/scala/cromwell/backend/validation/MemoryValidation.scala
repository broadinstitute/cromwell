package cromwell.backend.validation

import com.typesafe.config.Config
import squants.information.{Bytes, Information}
import wom.RuntimeAttributesKeys
import wom.values.{WomInteger, WomString}

import scala.util.{Failure, Success}

/**
  * Validates the "memory" runtime attribute as an Integer or String with format '8 GB', returning the value as a
  * `MemorySize`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * `configDefaultWdlValue` returns the value of the attribute as specified by the
  * reference.conf file, coerced into a WomValue.
  *
  * `optional` can be used to return the validated value as an `Option`,
  * wrapped in a `Some`, if present, or `None` if not found.
  *
  * `withDefaultMemory` can be used to create a memory validation that defaults to a particular memory size.
  */
object MemoryValidation {
  def instance(attributeName: String = RuntimeAttributesKeys.MemoryKey): RuntimeAttributesValidation[Information] =
    new MemoryValidation(attributeName)
  def optional(attributeName: String = RuntimeAttributesKeys.MemoryKey): OptionalRuntimeAttributesValidation[Information] =
    instance(attributeName).optional
  def configDefaultString(attributeName: String = RuntimeAttributesKeys.MemoryKey, config: Option[Config]): Option[String] =
    instance(attributeName).configDefaultValue(config)
  def withDefaultMemory(attributeName: String = RuntimeAttributesKeys.MemoryKey, memorySize: String): RuntimeAttributesValidation[Information] = {
    Information(memorySize) match {
      case Success(memory) => instance(attributeName).withDefault(WomInteger(memory.toBytes.toInt))
      case Failure(_) => instance(attributeName).withDefault(BadDefaultAttribute(WomString(memorySize.toString)))
    }
  }
}

class MemoryValidation(attributeName: String = RuntimeAttributesKeys.MemoryKey) extends InformationValidation(attributeName, Bytes)
