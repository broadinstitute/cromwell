package womtool.inputs

import cromwell.core.path.Path
import womtool.WomtoolMain.{SuccessfulTermination, Termination, UnsuccessfulTermination}
import womtool.input.WomGraphMaker
import wom.graph.{ExternalGraphInputNode, OptionalGraphInputNode, OptionalGraphInputNodeWithDefault, RequiredGraphInputNode}
import spray.json._
import spray.json.DefaultJsonProtocol._
import scala.util.{Failure, Success, Try}

object Inputs {
  def inputsJson(main: Path, showOptionals: Boolean): Termination = {

    WomGraphMaker.fromFiles(main, inputs = None) match {
      case Right(graph) =>
        Try(graph.externalInputNodes.toJson(inputNodeFormatter(showOptionals)).prettyPrint) match {
          case Success(json) => SuccessfulTermination(json + System.lineSeparator)
          case Failure(error) => UnsuccessfulTermination(error.getMessage)
        }
      case Left(errors) => UnsuccessfulTermination(errors.toList.mkString(System.lineSeparator))
    }
  }

  private def inputNodeFormatter(showOptionals: Boolean): JsonWriter[Set[ExternalGraphInputNode]] = set => {

    val valueMap: Seq[(String, JsValue)] = set.toList collect {
      case RequiredGraphInputNode(_, womType, nameInInputSet, _) => nameInInputSet -> JsString(womType.toDisplayString)
      case OptionalGraphInputNode(_, womOptionalType, nameInInputSet, _) if showOptionals => nameInInputSet -> JsString(womOptionalType.toDisplayString)
      case OptionalGraphInputNodeWithDefault(_, womType, default, nameInInputSet, _) if showOptionals => nameInInputSet -> JsString(s"${womType.toDisplayString} (default = ${default.sourceString})")
    }

    valueMap.toMap.toJson
  }
}
