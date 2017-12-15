package cromwell.webservice

import _root_.io.circe.yaml
import akka.util.ByteString
import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cats.instances.option._
import cats.syntax.apply._
import cats.syntax.functor._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.core._
import org.slf4j.LoggerFactory
import spray.json.{JsObject, JsValue}
import wdl.WorkflowJson
import wom.core._

import scala.util.{Failure, Success, Try}

final case class PartialWorkflowSources(workflowSource: WorkflowSource,
                                        workflowRoot: Option[String] = None,
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
  
  val WdlSourceKey = "wdlSource"
  val WorkflowRootKey = "workflowRoot"
  val WorkflowSourceKey = "workflowSource"
  val WorkflowTypeKey = "workflowType"
  val WorkflowTypeVersionKey = "workflowTypeVersion"
  val WorkflowInputsKey = "workflowInputs"
  val WorkflowInputsAuxPrefix = "workflowInputs_"
  val WorkflowOptionsKey = "workflowOptions"
  val labelsKey = "labels"
  val WdlDependenciesKey = "wdlDependencies"
  val WorkflowDependenciesKey = "workflowDependencies"

  val allKeys = List(WdlSourceKey, WorkflowRootKey, WorkflowSourceKey, WorkflowTypeKey, WorkflowTypeVersionKey, WorkflowInputsKey,
    WorkflowOptionsKey, labelsKey, WdlDependenciesKey, WorkflowDependenciesKey)

  val allPrefixes = List(WorkflowInputsAuxPrefix)

  def fromSubmitRoute(formData: Map[String, ByteString],
                      allowNoInputs: Boolean): Try[Seq[WorkflowSourceFilesCollection]] = {
    import cats.instances.list._
    import cats.syntax.apply._
    import cats.syntax.traverse._
    import cats.syntax.validated._

    val partialSources: ErrorOr[PartialWorkflowSources] = {
      def getStringValue(key: String) = formData.get(key).map(_.utf8String)
      def getArrayValue(key: String) = formData.get(key).map(_.toArray)

      // unrecognized keys
      val unrecognized: ErrorOr[Unit] = formData.keySet
        .filterNot(name => allKeys.contains(name) || allPrefixes.exists(name.startsWith))
        .toList
        .map(name => s"Unexpected body part name: $name")  match {
        case Nil => ().validNel
        case head :: tail => NonEmptyList.of(head, tail: _*).invalid
      }

      // workflow source
      val wdlSource = getStringValue(WdlSourceKey)
      val workflowSource = getStringValue(WorkflowSourceKey)

      def deprecationWarning(out: String, in: String)(actual: String): String = {
        if (actual == out) {
          val warning =
            Array(
              s"The '$out' parameter name has been deprecated in favor of '$in'.",
              s"Support for '$out' will be removed from future versions of Cromwell.",
              s"Please switch to using '$in' in future submissions.").mkString(" ")
          log.warn(warning)
          warning
        } else ""
      }

      val wdlSourceDeprecationWarning: String => String = deprecationWarning(out = WdlSourceKey, in = WorkflowSourceKey)
      val wdlSourceWarning = wdlSource.as(WdlSourceKey) map wdlSourceDeprecationWarning

      val workflowSourceFinal: ErrorOr[String] = (wdlSource, workflowSource) match {
        case (Some(source), None) => source.validNel
        case (None, Some(source)) => source.validNel
        case (Some(_), Some(_)) => s"$WdlSourceKey and $WorkflowSourceKey can't both be supplied".invalidNel
        case (None, None) => s"$WorkflowSourceKey needs to be supplied".invalidNel
      }

      // workflow inputs
      val workflowInputs: ErrorOr[Vector[WorkflowJson]] = getStringValue(WorkflowInputsKey) match {
        case Some(inputs) => workflowInputsValidation(inputs)
        case None => Vector.empty.validNel
      }
      val workflowInputsAux: ErrorOr[Map[Int, String]] = formData.toList.flatTraverse[ErrorOr, (Int, String)]({
        case (name, value) if name.startsWith(WorkflowInputsAuxPrefix) =>
          Try(name.stripPrefix(WorkflowInputsAuxPrefix).toInt).toErrorOr.map(index => List((index, value.utf8String)))
        case _ => List.empty.validNel
      }).map(_.toMap)

      // dependencies
      val wdlDependencies = getArrayValue(WdlDependenciesKey)
      val workflowDependencies = getArrayValue(WorkflowDependenciesKey)

      val wdlDependenciesDeprecationWarning: String => String = deprecationWarning(out = "wdlDependencies", in = "workflowDependencies")
      val wdlDependenciesWarning = wdlDependencies.as(WdlDependenciesKey) map wdlDependenciesDeprecationWarning

      val workflowDependenciesFinal: ErrorOr[Option[Array[Byte]]] = (wdlDependencies, workflowDependencies) match {
        case (Some(dep), None) => Option(dep).validNel
        case (None, Some(dep)) => Option(dep).validNel
        case (Some(_), Some(_)) => s"$WdlDependenciesKey and $WorkflowDependenciesKey can't both be supplied".invalidNel
        case (None, None) => None.validNel
      }

      (unrecognized, workflowSourceFinal, workflowInputs, workflowInputsAux, workflowDependenciesFinal) mapN {
        case (_, source, inputs, aux, dep) => PartialWorkflowSources(
          workflowSource = source,
          workflowRoot = getStringValue(WorkflowRootKey),
          workflowType = getStringValue(WorkflowTypeKey),
          workflowTypeVersion = getStringValue(WorkflowTypeVersionKey),
          workflowInputs = inputs,
          workflowInputsAux= aux,
          workflowOptions = getStringValue(WorkflowOptionsKey),
          customLabels = getStringValue(labelsKey),
          zippedImports = dep,
          warnings = wdlSourceWarning.toVector ++ wdlDependenciesWarning.toVector
        )
      }
    }

    partialSourcesToSourceCollections(partialSources, allowNoInputs) match {
      case Valid(source) => Success(source)
      case Invalid(errors) => Failure(new RuntimeException(s"Error(s): ${errors.toList.mkString("\n")}"))
    }
  }

  private def workflowInputsValidation(data: String): ErrorOr[Vector[WorkflowJson]] = {
    import _root_.io.circe.Printer
    import cats.syntax.validated._
    
    yaml.parser.parse(data) match {
      // If it's an array, treat each element as an individual input object, otherwise simply toString the whole thing
      case Right(json) => json.asArray.map(_.map(_.toString())).getOrElse(Vector(json.pretty(Printer.noSpaces))).validNel
      case Left(error) => s"Input file is not valid yaml nor json: ${error.getMessage}".invalidNel
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
        (validateInputs(partialSource),
          validateOptions(partialSource.workflowOptions), validateWorkflowType(partialSource),
          validateWorkflowTypeVersion(partialSource)) mapN {
          case (wfInputs, wfOptions, workflowType, workflowTypeVersion) =>
            wfInputs.map(inputsJson => WorkflowSourceFilesCollection(
              workflowSource = partialSource.workflowSource,
              workflowRoot = partialSource.workflowRoot,
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

