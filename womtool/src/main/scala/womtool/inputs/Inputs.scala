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
      case Right(graph) =>
        Try(graph.externalInputNodes.toJson(inputNodeWriter(showOptionals)).prettyPrint) match {
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

  private def womTypeToJson(womType: WomType, default: Option[WomExpression]): JsValue = womType match {
    case WomCompositeType(typeMap) => JsObject(typeMap.map {
      case (name, wt) => name -> womTypeToJson(wt, None)
    })
    case _ =>
      val defaultString = default.map(d => s"default = ${d.sourceString}").toList
      val optionalString = if (womType.isInstanceOf[WomOptionalType] || default.isDefined) List("optional") else List.empty

      val suffixStrings = optionalString ++ defaultString
      val suffixString = if (suffixStrings.nonEmpty) suffixStrings.mkString(" (", ", ", ")") else ""

      JsString(womType.toDisplayString + suffixString)
  }

}
