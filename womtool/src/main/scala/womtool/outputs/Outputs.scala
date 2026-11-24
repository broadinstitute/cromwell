package womtool.outputs

import cromwell.core.path.Path
import spray.json.{JsObject, JsString, JsValue}
import wom.graph.GraphOutputNode
import wom.types.{WomCompositeType, WomType}
import womtool.WomtoolMain.{SuccessfulTermination, Termination, UnsuccessfulTermination}
import womtool.input.WomGraphMaker

import scala.util.{Failure, Success, Try}

object Outputs {
  def outputsJson(main: Path): Termination =
    WomGraphMaker.fromFiles(main, inputs = None) match {
      case Right(graphWithImports) =>
        Try(outputNodesToJson(graphWithImports.graph.outputNodes).prettyPrint) match {
          case Success(json) => SuccessfulTermination(json + System.lineSeparator)
          case Failure(error) => UnsuccessfulTermination(error.getMessage)
        }
      case Left(errors) => UnsuccessfulTermination(errors.toList.mkString(System.lineSeparator))
    }

  private def outputNodesToJson(set: Set[GraphOutputNode]): JsValue =
    JsObject(set.toList.map(node => node.fullyQualifiedName -> womTypeToJson(node.womType)).toMap)

  private def womTypeToJson(womType: WomType): JsValue = womType match {
    case WomCompositeType(typeMap, _) =>
      JsObject(typeMap.map { case (name, wt) => name -> womTypeToJson(wt) })
    case _ => JsString(womType.stableName)
  }
}
