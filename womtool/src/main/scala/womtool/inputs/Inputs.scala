package womtool.inputs

import cromwell.core.path.Path
import womtool.WomtoolMain.{SuccessfulTermination, Termination, UnsuccessfulTermination}
import womtool.input.WomGraphMaker
import wom.graph.{ExternalGraphInputNode, OptionalGraphInputNode, OptionalGraphInputNodeWithDefault, RequiredGraphInputNode}
import spray.json._
import spray.json.DefaultJsonProtocol._
import wom.expression.WomExpression
import wom.types.{WomCompositeType, WomOptionalType, WomType}

import scala.util.{Failure, Success, Try}

object Inputs {
  def inputsJson(main: Path, showOptionals: Boolean): Termination = {

    WomGraphMaker.fromFiles(main, inputs = None) match {
      case Right(graphWithImports) =>
        Try(graphWithImports.graph.externalInputNodes.toJson(inputNodeWriter(showOptionals)).prettyPrint) match {
          case Success(json) => SuccessfulTermination(json + System.lineSeparator)
          case Failure(error) => UnsuccessfulTermination(error.getMessage)
        }
      case Left(errors) => UnsuccessfulTermination(errors.toList.mkString(System.lineSeparator))
    }
  }

  private def inputNodeWriter(showOptionals: Boolean): JsonWriter[Set[ExternalGraphInputNode]] = set => {

    val valueMap: Seq[(String, JsValue)] = set.toList collect {
      case RequiredGraphInputNode(_, womType, nameInInputSet, _) => nameInInputSet -> womTypeToJson(womType, None)
      case OptionalGraphInputNode(_, womOptionalType, nameInInputSet, _) if showOptionals => nameInInputSet -> womTypeToJson(womOptionalType, None)
      case OptionalGraphInputNodeWithDefault(_, womType, default, nameInInputSet, _) if showOptionals => nameInInputSet -> womTypeToJson(womType, Option(default))
    }

    valueMap.toMap.toJson
  }

  private def womTypeToJson(womType: WomType, default: Option[WomExpression]): JsValue = (womType, default) match {
    case (WomCompositeType(typeMap, _), _) => JsObject(
      typeMap.map { case (name, wt) => name -> womTypeToJson(wt, None) }
    )
    case (_, Some(d)) => JsString(s"${womType.stableName} (optional, default = ${d.sourceString})")
    case (_: WomOptionalType, _) => JsString(s"${womType.stableName} (optional)")
    case (_, _) => JsString(s"${womType.stableName}")
  }
}
