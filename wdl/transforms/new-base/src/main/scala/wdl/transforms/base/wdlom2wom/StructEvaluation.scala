package wdl.transforms.base.wdlom2wom

import cats.data.Validated.Valid
import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{StructElement, StructEntryElement}
import wom.types.{WomCompositeType, WomType}
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wdl.transforms.base.linking.typemakers._

object StructEvaluation {
  def convert(a: StructEvaluationInputs): ErrorOr[Map[String, WomType]] = {

    val validKnownAliases: Map[String, ErrorOr[WomType]] = a.knownTypeAliases map { case (key, value) => key -> value.validNel }
    val withNewStructs = a.structSections.foldLeft(validKnownAliases)(structFoldFunction)

    def unpackMapEntry(entry: (String, ErrorOr[WomType])): ErrorOr[(String, WomType)] = entry._2 map { entry._1 -> _ }

    withNewStructs.toList.traverse(unpackMapEntry).map(_.toMap)
  }

  private def structFoldFunction(current: Map[String, ErrorOr[WomType]], next: StructElement): Map[String, ErrorOr[WomType]] = {

    val currentValid = current collect {
      case (key, Valid(value)) => key -> value
    }

    def convertStructEntryElement(structEntryElement: StructEntryElement): ErrorOr[(String, WomType)] =
      structEntryElement.typeElement.determineWomType(currentValid).map((structEntryElement.identifier, _))

    val elementsValidation: ErrorOr[List[(String, WomType)]] = next.entries.toList.traverse(convertStructEntryElement)

    current + (next.name -> (elementsValidation map { elements => WomCompositeType(elements.toMap, Option(next.name)) } ))
  }

  case class StructEvaluationInputs(structSections: Seq[StructElement], knownTypeAliases: Map[String, WomType])
}

