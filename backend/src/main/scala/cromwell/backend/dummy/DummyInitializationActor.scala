package cromwell.backend.dummy

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
import cromwell.backend.validation.RuntimeAttributesValidation
import wom.expression.WomExpression
import wom.types.{WomStringType, WomType}
import wom.values.{WomString, WomValue}

class DummyInitializationActor(pipelinesParams: StandardInitializationActorParams)
  extends StandardInitializationActor(pipelinesParams) {

  override protected lazy val runtimeAttributeValidators: Map[String, Option[WomExpression] => Boolean] = Map("backend" -> { _ => true } )

  // Specific validator for "backend" to let me specify it in test cases (to avoid accidentally submitting the workflow to real backends!)
  val backendAttributeValidation: RuntimeAttributesValidation[String] = new RuntimeAttributesValidation[String] {
    override def key: String = "backend"

    override def coercion: Traversable[WomType] = Vector(WomStringType)

    override protected def validateValue: PartialFunction[WomValue, ErrorOr[String]] = {
      case WomString("Dummy") => "Dummy".validNel
      case other => s"Unexpected dummy backend value: $other".invalidNel
    }
  }

  override def runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder = super.runtimeAttributesBuilder.withValidation(backendAttributeValidation)
}
