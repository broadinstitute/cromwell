package cromwell.webservice

import _root_.io.circe.yaml
import akka.util.ByteString
import cats.Monoid
import cats.data.Validated.{Invalid, Valid}
import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.core._
import org.slf4j.LoggerFactory
import spray.json.{JsObject, JsValue}
import wdl.WorkflowJson
import wom.core._

import scala.util.Try

final case class PartialWorkflowSources(workflowSource: Option[WorkflowSource] = None,
                                        workflowType: Option[WorkflowType] = None,
                                        workflowTypeVersion: Option[WorkflowTypeVersion] = None,
                                        workflowInputs: Vector[WorkflowJson] = Vector.empty,
                                        workflowInputsAux: Map[Int, WorkflowJson] = Map.empty,
                                        workflowOptions: Option[WorkflowOptionsJson] = None,
                                        customLabels: Option[WorkflowJson] = None,
                                        zippedImports: Option[Array[Byte]] = None,
                                        warnings: Seq[String] = List.empty)

object PartialWorkflowSources {
  val log = LoggerFactory.getLogger(classOf[PartialWorkflowSources])

  implicit val partialWorkflowSourcesMonoid = new Monoid[PartialWorkflowSources] {
    override def empty: PartialWorkflowSources = PartialWorkflowSources()
    
    override def combine(x: PartialWorkflowSources, y: PartialWorkflowSources): PartialWorkflowSources = {
      x.copy(
        workflowSource = x.workflowSource.orElse(y.workflowSource),
        workflowType = x.workflowType.orElse(y.workflowType),
        workflowTypeVersion = x.workflowTypeVersion.orElse(y.workflowTypeVersion),
        workflowInputs = x.workflowInputs ++ y.workflowInputs,
        workflowInputsAux = x.workflowInputsAux ++ y.workflowInputsAux,
        workflowOptions = x.workflowOptions.orElse(y.workflowOptions),
        zippedImports = x.zippedImports.orElse(y.zippedImports),
        customLabels = x.customLabels.orElse(y.customLabels),
        warnings = x.warnings ++ y.warnings
      )
    }
  }
  
  def fromSubmitRoute(formData: Map[String, ByteString],
                      allowNoInputs: Boolean): Try[Seq[WorkflowSourceFilesCollection]] = {
    import cats.instances.list._
    import cats.syntax.foldable._
    import cats.syntax.validated._
    
    val partialSources: ErrorOr[PartialWorkflowSources] = {
      formData.toList.foldMap({
        case (name, data) =>
          name match {
            case "wdlSource" =>
              val warning = deprecationWarning(out = "wdlSource", in = "workflowSource")
              PartialWorkflowSources(workflowSource = Option(data.utf8String), warnings = Vector(warning)).validNel
            case "workflowSource" => 
              PartialWorkflowSources(workflowSource = Option(data.utf8String)).validNel
            case "workflowType" => 
              PartialWorkflowSources(workflowType = Option(data.utf8String)).validNel
            case "workflowTypeVersion" => PartialWorkflowSources(workflowTypeVersion = Option(data.utf8String)).validNel
            case "workflowInputs" =>
              workflowInputs(data.utf8String) map { inputs => PartialWorkflowSources(workflowInputs = inputs) }
            case _ if name.startsWith("workflowInputs_") =>
              Try(name.stripPrefix("workflowInputs_").toInt).toErrorOr map { index =>
                PartialWorkflowSources(workflowInputsAux = Map(index -> data.utf8String))
              }
            case "workflowOptions" => PartialWorkflowSources(workflowOptions = Option(data.utf8String)).validNel
            case "wdlDependencies" =>
              val warning = deprecationWarning(out = "wdlDependencies", in = "workflowDependencies")
              PartialWorkflowSources(zippedImports = Option(data.toArray), warnings = Vector(warning)).validNel
            case "workflowDependencies" => PartialWorkflowSources(zippedImports = Option(data.toArray)).validNel
            case "customLabels" => PartialWorkflowSources(customLabels = Option(data.utf8String)).validNel
            case _ => s"Unexpected body part name: $name".invalidNel
          }
      })
    }

    partialSourcesToSourceCollections(partialSources, allowNoInputs).toTry("Invalid submit request")
  }

  private def workflowInputs(data: String): ErrorOr[Vector[WorkflowJson]] = {
    import cats.syntax.validated._
    
    yaml.parser.parse(data) match {
      // If it's an array, treat each element as an individual input object, otherwise simply toString the whole thing
      case Right(json) => json.asArray.map(_.map(_.toString())).getOrElse(Vector(json.toString)).validNel
      case Left(error) => s"Input file is not valid json: ${error.getMessage}".invalidNel
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

