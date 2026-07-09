package womtool.inputs

import cromwell.core.path.Path
import spray.json.DefaultJsonProtocol._
import spray.json._
import wom.expression.WomExpression
import wom.graph._
import wom.types.{WomCompositeType, WomOptionalType, WomStringType, WomType}
import womtool.WomtoolMain.{SuccessfulTermination, Termination, UnsuccessfulTermination}
import womtool.input.WomGraphMaker

import scala.util.{Failure, Success, Try}

object Inputs {
  def inputsJson(main: Path, showOptionals: Boolean, showRuntimeOverrides: Boolean): Termination =
    WomGraphMaker.fromFiles(main, inputs = None) match {
      case Right(graphWithImports) =>
        Try(
          graphWithImports.graph.externalInputNodes
            .toJson(inputNodeWriter(showOptionals, showRuntimeOverrides))
            .prettyPrint
        ) match {
          case Success(json) => SuccessfulTermination(json + System.lineSeparator)
          case Failure(error) => UnsuccessfulTermination(error.getMessage)
        }
      case Left(errors) => UnsuccessfulTermination(errors.toList.mkString(System.lineSeparator))
    }

  private def inputNodeWriter(showOptionals: Boolean,
                              showRuntimeOverrides: Boolean
  ): JsonWriter[Set[ExternalGraphInputNode]] = set => {
    val valueMap: Seq[(String, JsValue)] =
      set.toList collect {
        case RequiredGraphInputNode(_, womType, nameInInputSet, _) => nameInInputSet -> womTypeToJson(womType, None)
        case OptionalGraphInputNode(_, womOptionalType, nameInInputSet, _) if showOptionals =>
          nameInInputSet -> womTypeToJson(womOptionalType, None)
        case OptionalGraphInputNodeWithDefault(_, womType, default, nameInInputSet, _) if showOptionals =>
          nameInInputSet -> womTypeToJson(womType, Option(default))
        // Special handling alert! Runtime override inputs are passed in individually as `$task.runtime.$attribute`
        // items, but they're combined internally to a single graph input node with a WomObjectType. Display them
        // in the `inputs` output with key `$task.runtime.*` to indicate that.
        case RuntimeOverrideGraphInputNode(_, nameInInputSet, _) if showRuntimeOverrides =>
          s"${nameInInputSet}.*" -> womTypeToJson(WomOptionalType(WomStringType), None)
      }

    valueMap.toMap.toJson
  }

  private def womTypeToJson(womType: WomType, default: Option[WomExpression]): JsValue = (womType, default) match {
    case (WomCompositeType(typeMap, _), _) =>
      JsObject(
        typeMap.map { case (name, wt) => name -> womTypeToJson(wt, None) }
      )
    case (_, Some(d)) => JsString(s"${womType.stableName} (optional, default = ${d.sourceString})")
    case (_: WomOptionalType, _) => JsString(s"${womType.stableName} (optional)")
    case (_, _) => JsString(s"${womType.stableName}")
  }
}
