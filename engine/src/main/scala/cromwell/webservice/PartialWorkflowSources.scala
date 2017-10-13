package cromwell.webservice

import akka.util.ByteString
import cats.data.Validated.{Invalid, Valid}
import cats.syntax.apply._
import cats.syntax.validated._
import cromwell.core._
import lenthall.validation.ErrorOr.ErrorOr
import lenthall.validation.Validation._
import org.slf4j.LoggerFactory
import spray.json.{JsObject, JsValue}
import wdl.WorkflowJson
import wom.core._

import scala.util.Try

final case class PartialWorkflowSources(workflowSource: Option[WorkflowSource],
                                        workflowType: Option[WorkflowType],
                                        workflowTypeVersion: Option[WorkflowTypeVersion],
                                        workflowInputs: Vector[WorkflowJson],
                                        workflowInputsAux: Map[Int, WorkflowJson],
                                        workflowOptions: Option[WorkflowOptionsJson],
                                        customLabels: Option[WorkflowJson],
                                        zippedImports: Option[Array[Byte]],
                                        warnings: Seq[String])

object PartialWorkflowSources {
  val log = LoggerFactory.getLogger(classOf[PartialWorkflowSources])

  def empty = PartialWorkflowSources(
    workflowSource = None,
    workflowType = None,
    workflowTypeVersion = None,
    workflowInputs = Vector.empty,
    workflowInputsAux = Map.empty,
    workflowOptions = None,
    customLabels = None,
    zippedImports = None,
    warnings = Vector.empty
  )

  def fromSubmitRoute(formData: Map[String, ByteString],
                      allowNoInputs: Boolean): Try[Seq[WorkflowSourceFilesCollection]] = {
    val partialSources: Try[PartialWorkflowSources] = Try {
      formData.foldLeft(PartialWorkflowSources.empty) {
        (partialSources: PartialWorkflowSources, kv: (String, ByteString)) =>
          val (name, data) = kv

          name match {
            case "wdlSource" =>
              val warning = deprecationWarning(out = "wdlSource", in = "workflowSource")
              val warnings = warning +: partialSources.warnings
              partialSources.copy(workflowSource = Option(data.utf8String), warnings = warnings)
            case "workflowSource" => partialSources.copy(workflowSource = Option(data.utf8String))
            case "workflowType" => partialSources.copy(workflowType = Option(data.utf8String))
            case "workflowTypeVersion" => partialSources.copy(workflowTypeVersion = Option(data.utf8String))
            case "workflowInputs" => partialSources.copy(workflowInputs = workflowInputs(data.utf8String))
            case _ if name.startsWith("workflowInputs_") =>
              val index = name.stripPrefix("workflowInputs_").toInt
              partialSources.copy(workflowInputsAux = partialSources.workflowInputsAux + (index -> data.utf8String))
            case "workflowOptions" => partialSources.copy(workflowOptions = Option(data.utf8String))
            case "wdlDependencies" =>
              val warning = deprecationWarning(out = "wdlDependencies", in = "workflowDependencies")
              val warnings = warning +: partialSources.warnings
              partialSources.copy(zippedImports = Option(data.toArray), warnings = warnings)
            case "workflowDependencies" => partialSources.copy(zippedImports = Option(data.toArray))
            case "customLabels" => partialSources.copy(customLabels = Option(data.utf8String))
            case _ => throw new IllegalArgumentException(s"Unexpected body part name: $name")
          }
      }
    }

    partialSourcesToSourceCollections(partialSources.toErrorOr, allowNoInputs).toTry
  }

  private def workflowInputs(data: String): Vector[WorkflowJson] = {
    import spray.json._
    data.parseJson match {
      case JsArray(Seq(x, xs@_*)) => (Vector(x) ++ xs).map(_.compactPrint)
      case JsArray(_) => Vector.empty
      case v: JsValue => Vector(v.compactPrint)
    }
  }

  private def partialSourcesToSourceCollections(partialSources: ErrorOr[PartialWorkflowSources],
                                                allowNoInputs: Boolean): ErrorOr[Seq[WorkflowSourceFilesCollection]] = {
    def validateInputs(pws: PartialWorkflowSources): ErrorOr[Seq[WorkflowJson]] =
      (pws.workflowInputs.isEmpty, allowNoInputs) match {
        case (true, true) => Vector("{}").validNel
        case (true, false) => "No inputs were provided".invalidNel
        case _ =>
          val sortedInputAuxes = pws.workflowInputsAux.toSeq.sortBy { case (index, _) => index } map { case(_, inputJson) => Option(inputJson) }
          (pws.workflowInputs map { workflowInputSet: WorkflowJson => mergeMaps(Seq(Option(workflowInputSet)) ++ sortedInputAuxes).toString }).validNel
      }

    def validateOptions(options: Option[WorkflowOptionsJson]): ErrorOr[WorkflowOptions] =
      WorkflowOptions.fromJsonString(options.getOrElse("{}")).toErrorOr leftMap { _ map { i => s"Invalid workflow options provided: $i" } }

    def validateWorkflowSource(partialSource: PartialWorkflowSources): ErrorOr[WorkflowJson] = partialSource.workflowSource match {
      case Some(src) => src.validNel
      case _ => s"Incomplete workflow submission: $partialSource".invalidNel
    }

    def validateWorkflowType(partialSource: PartialWorkflowSources): ErrorOr[Option[WorkflowType]] = {
      partialSource.workflowType match {
        case Some(_) => partialSource.workflowType.validNel
        case None => WorkflowOptions.defaultWorkflowType.validNel
      }
    }

    def validateWorkflowTypeVersion(partialSource: PartialWorkflowSources): ErrorOr[Option[WorkflowTypeVersion]] = {
      partialSource.workflowTypeVersion match {
        case Some(src) => Option(src).validNel
        case None => WorkflowOptions.defaultWorkflowTypeVersion.validNel
      }
    }

    partialSources match {
      case Valid(partialSource) =>
        (validateWorkflowSource(partialSource), validateInputs(partialSource),
          validateOptions(partialSource.workflowOptions), validateWorkflowType(partialSource),
          validateWorkflowTypeVersion(partialSource)) mapN {
          case (wfSource, wfInputs, wfOptions, workflowType, workflowTypeVersion) =>
            wfInputs.map(inputsJson => WorkflowSourceFilesCollection(
              workflowSource = wfSource,
              workflowType = workflowType,
              workflowTypeVersion = workflowTypeVersion,
              inputsJson = inputsJson,
              workflowOptionsJson = wfOptions.asPrettyJson,
              labelsJson = partialSource.customLabels.getOrElse("{}"),
              importsFile = partialSource.zippedImports,
              warnings = partialSource.warnings))
        }
      case Invalid(err) => err.invalid
    }
  }

  private def deprecationWarning(out: String, in: String): String = {
    val warning =
      s"""
         |The '$out' parameter name has been deprecated in favor of '$in'.
         |Support for '$out' will be removed from future versions of Cromwell.
         |Please switch to using '$in' in future submissions.
         """.stripMargin
    log.warn(warning)
    warning
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

