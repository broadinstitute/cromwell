package womtool.parametermeta

import common.Checked
import cromwell.core.path.Path
import spray.json.{JsNull, JsObject, JsonWriter, enrichAny}
import spray.json.DefaultJsonProtocol._
import wom.callable.{MetaValueElement, TaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle
import womtool.WomtoolMain.{SuccessfulTermination, Termination, UnsuccessfulTermination}
import womtool.input.WomGraphMaker

object ParameterMeta {
  def parameterMetaInfoJson(main: Path): Termination = {
    WomGraphMaker.getBundle(main) match {
      case Right(b) =>
        getParameterMeta(b) match {
          case Right(inputs) =>
            SuccessfulTermination(inputs.toJson(parameterMetaInfoJsonWriter()).prettyPrint)
          case Left(errors) =>
            UnsuccessfulTermination(errors.toList.mkString(System.lineSeparator))
        }
      case Left(errors) =>
        UnsuccessfulTermination(errors.toList.mkString(System.lineSeparator))
    }
  }

  private final case class ParameterMetaInfo(inputs: Map[String, Option[MetaValueElement]],
                                             outputs: Map[String, Option[MetaValueElement]])

  private def parameterMetaInfoJsonWriter(): JsonWriter[ParameterMetaInfo] = info => {
    val inputs = info.inputs.toJson(parameterMetaJsonWriter())
    val outputs = info.outputs.toJson(parameterMetaJsonWriter())
    JsObject(
      "inputs" -> inputs,
      "outputs" -> outputs
    )
  }

  private def parameterMetaJsonWriter(): JsonWriter[Map[String, Option[MetaValueElement]]] = m => {
    m.collect {
      case (name, Some(meta)) => name -> meta.toJson(metaValueElementJsonWriter())
    }.toJson
  }

  private def metaValueElementJsonWriter(): JsonWriter[MetaValueElement] = {
    case MetaValueElement.MetaValueElementNull => JsNull
    case MetaValueElement.MetaValueElementBoolean(value) => value.toJson
    case MetaValueElement.MetaValueElementFloat(value) => value.toJson
    case MetaValueElement.MetaValueElementInteger(value) => value.toJson
    case MetaValueElement.MetaValueElementString(value) => value.toJson
    case MetaValueElement.MetaValueElementObject(value) =>
      value.view.mapValues(m => m.toJson(metaValueElementJsonWriter())).toMap.toJson
    case MetaValueElement.MetaValueElementArray(value) =>
      value.map(m => m.toJson(metaValueElementJsonWriter())).toJson
  }

  private def getParameterMeta(b: WomBundle): Checked[ParameterMetaInfo] = {
    b.toExecutableCallable.map(primaryCallable => {
      // inputs the same as womtool inputs
      val inputs = primaryCallable.graph.externalInputNodes.map(inputNode => {
        val parameterFullName = inputNode.nameInInputSet
        parameterFullName -> findParameterMetaValue(b, parameterFullName)
      }).toMap
      // outputs the same as womtool outputs
      val outputs = primaryCallable.graph.outputNodes.map(outputNode => {
        val parameterFullName = outputNode.fullyQualifiedName
        parameterFullName -> findParameterMetaValue(b, parameterFullName)
      }).toMap
      ParameterMetaInfo(inputs, outputs)
    })
  }

  private def findParameterMetaValue(b: WomBundle, parameterFullName: String): Option[MetaValueElement] = {
    val parts = parameterFullName.split("\\.")
    val parameterName = parts(parts.length - 1)
    val callableName: Option[String] = if (parts.length >= 2) Some(parts(parts.length - 2)) else None
    callableName.flatMap(name => {
      b.allCallables.get(name) match {
        case Some(w: WorkflowDefinition) => w.parameterMeta.get(parameterName)
        case Some(t: TaskDefinition) => t.parameterMeta.get(parameterName)
        case _ => None
      }
    })
  }

}