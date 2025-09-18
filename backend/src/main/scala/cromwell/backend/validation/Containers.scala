package cromwell.backend.validation

import cats.implicits.catsSyntaxValidatedId
import common.validation.ErrorOr.ErrorOr
import wom.RuntimeAttributesKeys
import wom.types.{WomStringType, WomType}
import wom.values.{WomArray, WomString, WomValue}

// Wrapper type for the list of containers that can be provided as 'docker' or 'container' runtime attributes.
// In WDL, this value can be either a single string or an array of strings, but in the backend we always
// want to deal with it as a list of strings.
//
// Previous to WDL 1.1, only 'docker' was supported. From WDL 1.1 onwards, they are aliases of each other, with
// `container` being preferred and `docker` deprecated. Only one of they two may be provided in runtime attrs.
// Note that we strip `container` out of pre-1.1 WDL files during parsing, so at this stage we only see `docker`
// in those cases.
case class Containers(values: List[String])

object Containers {
  val validWdlTypes: Set[wom.types.WomType] =
    Set(wom.types.WomStringType, wom.types.WomArrayType(wom.types.WomStringType))

  val runtimeAttrKeys = Set(RuntimeAttributesKeys.DockerKey, RuntimeAttributesKeys.ContainerKey)

  def apply(value: String): Containers = Containers(List(value))

  def extractContainer(validatedRuntimeAttributes: ValidatedRuntimeAttributes): String =
    extractContainerOption(validatedRuntimeAttributes).getOrElse {
      throw new RuntimeException("No container image found in either 'container' or 'docker' runtime attributes.")
    }

  def extractContainerOption(validatedRuntimeAttributes: ValidatedRuntimeAttributes): Option[String] = {
    val dockerContainer = RuntimeAttributesValidation
      .extractOption(DockerValidation.instance, validatedRuntimeAttributes)
      .flatMap(_.values.headOption)
    val containerContainer = RuntimeAttributesValidation
      .extractOption(ContainerValidation.instance, validatedRuntimeAttributes)
      .flatMap(_.values.headOption)

    containerContainer.orElse(dockerContainer)
  }
}

/**
 * Trait to handle validation of both 'docker' and 'container' runtime attributes, which are mutually exclusive
 * ways of specifying the container image to use for a task.
 */
trait ContainersValidation extends RuntimeAttributesValidation[Containers] {
  override def coercion: Set[WomType] = Containers.validWdlTypes

  override def usedInCallCaching: Boolean = true

  override protected def missingValueMessage: String = s"Can't find an attribute value for key ${key}"

  override protected def invalidValueMessage(value: WomValue): String = super.missingValueMessage

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Containers]] = {
    case WomString(value) => value.validNel.map(v => Containers(v))
    case WomArray(womType, values) if womType.memberType == WomStringType =>
      Containers(values.map(_.valueString).toList).validNel
  }

  override def validate(values: Map[String, WomValue]): ErrorOr[Containers] =
    if (Containers.runtimeAttrKeys.count(v => values.contains(v)) > 1) {
      s"Must provide only one of the following runtime attributes: ${Containers.runtimeAttrKeys.mkString(", ")}".invalidNel
    } else super.validate(values)
}
