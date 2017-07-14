package cromwell.webservice

import akka.util.ByteString
import cromwell.core.{WorkflowOptions, WorkflowOptionsJson, WorkflowSourceFilesCollection}
import wdl4s.wdl.{WorkflowJson, WorkflowSource}
import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import cats.syntax.cartesian._
import lenthall.validation.ErrorOr.ErrorOr
import cromwell.core._
import org.slf4j.LoggerFactory
import spray.json.{JsObject, JsValue}

import scala.util.Try

final case class PartialWorkflowSources(workflowSource: Option[WorkflowSource],
                                        workflowType: Option[WorkflowType],
                                        workflowTypeVersion: Option[WorkflowTypeVersion],
                                        workflowInputs: Vector[WorkflowJson],
                                        workflowInputsAux: Map[Int, WorkflowJson],
                                        workflowOptions: Option[WorkflowOptionsJson],
                                        customLabels: Option[WorkflowJson],
                                        zippedImports: Option[Array[Byte]])

object PartialWorkflowSources {
  val log = LoggerFactory.getLogger(classOf[PartialWorkflowSources])

  def empty = PartialWorkflowSources(
    workflowSource = None,
    // TODO do not hardcode, especially not out here at the boundary layer good gravy
    workflowType = Option("WDL"),
    workflowTypeVersion = None,
    workflowInputs = Vector.empty,
    workflowInputsAux = Map.empty,
    workflowOptions = None,
    customLabels = None,
    zippedImports = None
  )

  def fromSubmitRoute(formData: Map[String, ByteString], allowNoInputs: Boolean): Try[Seq[WorkflowSourceFilesCollection]] = {
    val partialSources = Try(formData.foldLeft(PartialWorkflowSources.empty) { (partialSources: PartialWorkflowSources, kv: (String, ByteString)) =>
      val name = kv._1
      val data = kv._2

      if (name == "WorkflowSource" || name == "workflowSource") {
        if (name == "WorkflowSource") deprecationWarning(out = "WorkflowSource", in = "workflowSource")
        partialSources.copy(workflowSource = Option(data.utf8String))
      } else if (name == "workflowType") {
        partialSources.copy(workflowType = Option(data.utf8String))
      } else if (name == "workflowTypeVersion") {
        partialSources.copy(workflowTypeVersion = Option(data.utf8String))
      } else if (name == "workflowInputs") {
        partialSources.copy(workflowInputs = workflowInputs(data.utf8String))
      } else if (name.startsWith("workflowInputs_")) {
        val index = name.stripPrefix("workflowInputs_").toInt
        partialSources.copy(workflowInputsAux = partialSources.workflowInputsAux + (index -> data.utf8String))
      } else if (name == "workflowOptions") {
        partialSources.copy(workflowOptions = Option(data.utf8String))
      } else if (name == "wdlDependencies" || name == "workflowDependencies") {
        if (name == "wdlDependencies") deprecationWarning(out = "wdlDependencies", in = "workflowDependencies")
        partialSources.copy(zippedImports = Option(data.toArray))
      } else if (name == "customLabels") {
        partialSources.copy(customLabels = Option(data.utf8String))
      } else {
        throw new IllegalArgumentException(s"Unexpected body part name: $name")
      }
    })

    partialSourcesToSourceCollections(partialSources.tryToErrorOr, allowNoInputs).errorOrToTry
  }

  private def workflowInputs(data: String): Vector[WorkflowJson] = {
    import spray.json._
    data.parseJson match {
      case JsArray(Seq(x, xs@_*)) => (Vector(x) ++ xs).map(_.compactPrint)
      case JsArray(_) => Vector.empty
      case v: JsValue => Vector(v.compactPrint)
    }
  }

  private def partialSourcesToSourceCollections(partialSources: ErrorOr[PartialWorkflowSources], allowNoInputs: Boolean): ErrorOr[Seq[WorkflowSourceFilesCollection]] = {
    def validateInputs(pws: PartialWorkflowSources): ErrorOr[Seq[WorkflowJson]] =
      (pws.workflowInputs.isEmpty, allowNoInputs) match {
        case (true, true) => Vector("{}").validNel
        case (true, false) => "No inputs were provided".invalidNel
        case _ =>
          val sortedInputAuxes = pws.workflowInputsAux.toSeq.sortBy { case (index, _) => index } map { case(_, inputJson) => Option(inputJson) }
          (pws.workflowInputs map { workflowInputSet: WorkflowJson => mergeMaps(Seq(Option(workflowInputSet)) ++ sortedInputAuxes).toString }).validNel
      }

    def validateOptions(options: Option[WorkflowOptionsJson]): ErrorOr[WorkflowOptions] =
      WorkflowOptions.fromJsonString(options.getOrElse("{}")).tryToErrorOr leftMap { _ map { i => s"Invalid workflow options provided: $i" } }

    def validateWorkflowSource(partialSource: PartialWorkflowSources): ErrorOr[WorkflowJson] = partialSource.workflowSource match {
      case Some(src) => src.validNel
      case _ => s"Incomplete workflow submission: $partialSource".invalidNel
    }

    partialSources match {
      case Valid(partialSource) =>
        (validateWorkflowSource(partialSource) |@| validateInputs(partialSource) |@| validateOptions(partialSource.workflowOptions)) map {
          case (wfSource, wfInputs, wfOptions) =>
            wfInputs.map(inputsJson => WorkflowSourceFilesCollection(
              workflowSource = wfSource,
              workflowType = partialSource.workflowType,
              workflowTypeVersion = partialSource.workflowTypeVersion,
              inputsJson = inputsJson,
              workflowOptionsJson = wfOptions.asPrettyJson,
              labelsJson = partialSource.customLabels.getOrElse("{}"),
              importsFile = partialSource.zippedImports))        }
      case Invalid(err) => err.invalid
    }
  }

  private def deprecationWarning(out: String, in: String): Unit = {
    val warning =
      s"""
         |The '$out' parameter name has been deprecated in favor of '$in'.
         |Support for '$out' will be removed from future versions of Cromwell.
         |Please switch to using '$in' in future submissions.
         """.stripMargin
    log.warn(warning)
  }

  def mergeMaps(allInputs: Seq[Option[String]]): JsObject = {
    val convertToMap = allInputs.map(x => toMap(x))
    JsObject(convertToMap reduce (_ ++ _))
  }

  private def toMap(someInput: Option[String]): Map[String, JsValue] = {
    import spray.json._
    someInput match {
      case Some(inputs: String) => inputs.parseJson match {
        case JsObject(inputMap) => inputMap
        case _ =>
          throw new RuntimeException(s"Submitted inputs couldn't be processed, please check for syntactical errors")
      }
      case None => Map.empty
    }
  }
}

