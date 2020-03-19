package womtool.outputs

import cromwell.core.path.Path
import spray.json.DefaultJsonProtocol._
import spray.json._
import wom.graph.GraphOutputNode
import wom.types.{WomCompositeType, WomType}
import womtool.WomtoolMain.{SuccessfulTermination, Termination, UnsuccessfulTermination}
import womtool.input.WomGraphMaker

import scala.util.{Failure, Success, Try}

object Outputs {
  def outputsJson(main: Path): Termination = {

    WomGraphMaker.fromFiles(main, inputs = None) match {
      case Right(graphWithImports) =>
        Try(graphWithImports.graph.outputNodes.toJson(outputNodeWriter()).prettyPrint) match {
          case Success(json) => SuccessfulTermination(json + System.lineSeparator)
          case Failure(error) => UnsuccessfulTermination(error.getMessage)
        }
      case Left(errors) => UnsuccessfulTermination(errors.toList.mkString(System.lineSeparator))
    }
  }

  private def outputNodeWriter(): JsonWriter[Set[GraphOutputNode]] = set => {

    val valueMap: Seq[(String, JsValue)] = set.toList collect {
      case node: GraphOutputNode => node.fullyQualifiedName -> womTypeToJson(node.womType)
    }

    valueMap.toMap.toJson
  }

  private def womTypeToJson(womType: WomType): JsValue = womType match {
    case WomCompositeType(typeMap, _) => JsObject(
      typeMap.map { case (name, wt) => name -> womTypeToJson(wt) }
    )
    case _ => JsString(s"${womType.stableName}")
  }
}
