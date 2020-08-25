package womtool.outputs

import io.circe.syntax._
import cromwell.core.path.Path
import io.circe.{Encoder, Json, JsonObject, Printer}
import wom.graph.GraphOutputNode
import wom.types.{WomCompositeType, WomType}
import womtool.WomtoolMain.{SuccessfulTermination, Termination, UnsuccessfulTermination}
import womtool.input.WomGraphMaker

import scala.util.{Failure, Success, Try}

object Outputs {
  def outputsJson(main: Path): Termination = {

    WomGraphMaker.fromFiles(main, inputs = None) match {
      case Right(graphWithImports) =>
        Try(graphWithImports.graph.outputNodes.asJson.printWith(sprayLikePrettyPrinter)) match {
          case Success(json) => SuccessfulTermination(json + System.lineSeparator)
          case Failure(error) => UnsuccessfulTermination(error.getMessage)
        }
      case Left(errors) => UnsuccessfulTermination(errors.toList.mkString(System.lineSeparator))
    }
  }

  private implicit val e: Encoder[Set[GraphOutputNode]] = (set: Set[GraphOutputNode]) => {
    val valueMap = set.toList.map {
      node: GraphOutputNode => node.fullyQualifiedName -> womTypeToJson(node.womType)
    }

    valueMap.toMap.asJson
  }

  private def womTypeToJson(womType: WomType): Json = womType match {
    case WomCompositeType(typeMap, _) => JsonObject.fromMap(
      typeMap.map { case (name, wt) => name -> womTypeToJson(wt) }
    ).asJson
    case _ => womType.stableName.asJson
  }

  private val sprayLikePrettyPrinter = Printer.spaces2.copy(dropNullValues = true, colonLeft = "")
}
