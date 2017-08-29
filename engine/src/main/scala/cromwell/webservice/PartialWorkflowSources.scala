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
                      allowNoInputs: Boolean,
                      version: Int): Try[Seq[WorkflowSourceFilesCollection]] = {
    val partialSources: Try[PartialWorkflowSources] = Try {
      formData.foldLeft(PartialWorkflowSources.empty) {
        (partialSources: PartialWorkflowSources, kv: (String, ByteString)) =>
          val (name, data) = kv

          (name, version) match {
            case ("wdlSource", 1) =>
              val warning = deprecationWarning(out = "wdlSource", in = "workflowSource")
              val warnings = warning +: partialSources.warnings
              partialSources.copy(workflowSource = Option(data.utf8String), warnings = warnings)
            case ("workflowSource", _) => partialSources.copy(workflowSource = Option(data.utf8String))
            case ("workflowType", _) => partialSources.copy(workflowType = Option(data.utf8String))
            case ("workflowTypeVersion", _) => partialSources.copy(workflowTypeVersion = Option(data.utf8String))
            case ("workflowInputs", _) => partialSources.copy(workflowInputs = workflowInputs(data.utf8String))
            case (_, _) if name.startsWith("workflowInputs_") =>
              val index = name.stripPrefix("workflowInputs_").toInt
              partialSources.copy(workflowInputsAux = partialSources.workflowInputsAux + (index -> data.utf8String))
            case ("workflowOptions", _) => partialSources.copy(workflowOptions = Option(data.utf8String))
            case ("wdlDependencies", 1) =>
              val warning = deprecationWarning(out = "wdlDependencies", in = "workflowDependencies")
              val warnings = warning +: partialSources.warnings
              partialSources.copy(zippedImports = Option(data.toArray), warnings = warnings)
            case ("workflowDependencies", _) => partialSources.copy(zippedImports = Option(data.toArray))
            case ("customLabels", _) => partialSources.copy(customLabels = Option(data.utf8String))
            case _ => throw new IllegalArgumentException(s"Unexpected body part name: $name")
          }
      }
    }

    partialSourcesToSourceCollections(partialSources.tryToErrorOr, allowNoInputs, version).errorOrToTry
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
                                                allowNoInputs: Boolean,
                                                version: Int): ErrorOr[Seq[WorkflowSourceFilesCollection]] = {
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

    def validateWorkflowType(partialSource: PartialWorkflowSources, version: Int): ErrorOr[Option[WorkflowType]] = {
      (partialSource.workflowType, version) match {
        case (Some(src), _) => Option(src).validNel
        case (None, 1) => Option("WDL").validNel
        case (None, _) => s"Workflow type is mandatory for v$version".invalidNel
      }
    }

    def validateWorkflowTypeVersion(partialSource: PartialWorkflowSources,
                                    version: Int): ErrorOr[Option[WorkflowTypeVersion]] = {
      (partialSource.workflowTypeVersion, version) match {
        case (Some(src), _) => Option(src).validNel
        case (None, 1) => None.validNel
        case (None, _) => s"Workflow type version is mandatory for v$version".invalidNel
      }
    }

    partialSources match {
      case Valid(partialSource) =>
        (validateWorkflowSource(partialSource) |@|
          validateInputs(partialSource) |@|
          validateOptions(partialSource.workflowOptions) |@|
          validateWorkflowType(partialSource, version) |@|
          validateWorkflowTypeVersion(partialSource, version)) map {
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

