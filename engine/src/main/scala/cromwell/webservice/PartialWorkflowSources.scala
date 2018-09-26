package cromwell.webservice

import java.net.URL

import _root_.io.circe.yaml
import akka.util.ByteString
import cats.data.NonEmptyList
import cats.data.Validated._
import cats.instances.list._
import cats.instances.option._
import cats.syntax.apply._
import cats.syntax.functor._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.core._
import cromwell.core.labels.Label
import org.slf4j.LoggerFactory
import spray.json.{JsObject, JsValue, _}
import wdl.draft2.model.WorkflowJson
import wom.core._

import scala.util.{Failure, Success, Try}

final case class PartialWorkflowSources(workflowSource: Option[WorkflowSource] = None,
                                        workflowUrl: Option[WorkflowUrl] = None,
                                        workflowRoot: Option[String] = None,
                                        workflowType: Option[WorkflowType] = None,
                                        workflowTypeVersion: Option[WorkflowTypeVersion] = None,
                                        workflowInputs: Vector[WorkflowJson] = Vector.empty,
                                        workflowInputsAux: Map[Int, WorkflowJson] = Map.empty,
                                        workflowOptions: Option[WorkflowOptionsJson] = None,
                                        customLabels: Option[WorkflowJson] = None,
                                        zippedImports: Option[Array[Byte]] = None,
                                        warnings: Seq[String] = List.empty,
                                        workflowOnHold: Boolean)

object PartialWorkflowSources {
  val log = LoggerFactory.getLogger(classOf[PartialWorkflowSources])

  val WdlSourceKey = "wdlSource"
  val WorkflowRootKey = "workflowRoot"
  val WorkflowSourceKey = "workflowSource"
  val WorkflowUrlKey = "workflowUrl"
  val WorkflowTypeKey = "workflowType"
  val WorkflowTypeVersionKey = "workflowTypeVersion"
  val WorkflowInputsKey = "workflowInputs"
  val WorkflowInputsAuxPrefix = "workflowInputs_"
  val WorkflowOptionsKey = "workflowOptions"
  val labelsKey = "labels"
  val WdlDependenciesKey = "wdlDependencies"
  val WorkflowDependenciesKey = "workflowDependencies"
  val workflowOnHoldKey = "workflowOnHold"

  val allKeys = List(WdlSourceKey, WorkflowUrlKey, WorkflowRootKey, WorkflowSourceKey, WorkflowTypeKey, WorkflowTypeVersionKey, WorkflowInputsKey,
    WorkflowOptionsKey, labelsKey, WdlDependenciesKey, WorkflowDependenciesKey, workflowOnHoldKey)

  val allPrefixes = List(WorkflowInputsAuxPrefix)

  val MaxWorkflowUrlLength = 2000

  def fromSubmitRoute(formData: Map[String, ByteString],
                      allowNoInputs: Boolean): Try[Seq[WorkflowSourceFilesCollection]] = {
    import cats.instances.list._
    import cats.syntax.apply._
    import cats.syntax.traverse._
    import cats.syntax.validated._

    val partialSources: ErrorOr[PartialWorkflowSources] = {
      def getStringValue(key: String) = formData.get(key).map(_.utf8String)
      def getBooleanValue(key: String) = getStringValue(key).map(b => Try(b.toBoolean).toErrorOr)
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
      val workflowUrl = getStringValue(WorkflowUrlKey)

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

      val workflowSourceFinal: ErrorOr[Option[String]] = (wdlSource, workflowSource, workflowUrl) match {
        case (Some(source), None, None) => Some(source).validNel
        case (None, Some(source), None) => Some(source).validNel
        case (None, None, Some(_)) => None.validNel
        case (Some(_), Some(_), None) => s"$WdlSourceKey and $WorkflowSourceKey can't both be supplied".invalidNel
        case (None, Some(_), Some(_)) => s"$WorkflowSourceKey and $WorkflowUrlKey can't both be supplied".invalidNel
        case (Some(_), None, Some(_)) => s"$WdlSourceKey and $WorkflowUrlKey can't both be supplied".invalidNel
        case (Some(_), Some(_), Some(_)) => s"$WdlSourceKey, $WorkflowSourceKey and $WorkflowUrlKey all 3 can't be supplied".invalidNel
        case (None, None, None) => s"$WorkflowSourceKey or $WorkflowUrlKey needs to be supplied".invalidNel
      }

      // workflow inputs
      val workflowInputs: ErrorOr[Vector[WorkflowJson]] = getStringValue(WorkflowInputsKey) match {
        case Some(inputs) => workflowInputsValidation(inputs)
        case None => Vector.empty.validNel
      }
      val workflowInputsAux: ErrorOr[Map[Int, String]] = formData.toList.flatTraverse({
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

      val onHold: ErrorOr[Boolean] = getBooleanValue(workflowOnHoldKey).getOrElse(false.validNel)

      (unrecognized, workflowSourceFinal, workflowInputs, workflowInputsAux, workflowDependenciesFinal, onHold) mapN {
        case (_, source, inputs, aux, dep, onHold) => PartialWorkflowSources(
          workflowSource = source,
          workflowUrl = workflowUrl,
          workflowRoot = getStringValue(WorkflowRootKey),
          workflowType = getStringValue(WorkflowTypeKey),
          workflowTypeVersion = getStringValue(WorkflowTypeVersionKey),
          workflowInputs = inputs,
          workflowInputsAux= aux,
          workflowOptions = getStringValue(WorkflowOptionsKey),
          customLabels = getStringValue(labelsKey),
          zippedImports = dep,
          warnings = wdlSourceWarning.toVector ++ wdlDependenciesWarning.toVector,
          workflowOnHold = onHold
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

    def validateLabels(labels: WorkflowJson) : ErrorOr[WorkflowJson] = {

      def validateKeyValuePair(key: String, value: String): ErrorOr[Unit] = (Label.validateLabelKey(key), Label.validateLabelValue(value)).tupled.void

      def validateLabelRestrictions(inputs: Map[String, JsValue]): ErrorOr[Unit] = {
        inputs.toList.traverse[ErrorOr, Unit]({
          case (key, JsString(s)) => validateKeyValuePair(key, s)
          case (key, other) => s"Invalid label $key: $other : Labels must be strings. ${Label.LabelExpectationsMessage}".invalidNel
        }).void
      }


      Try(labels.parseJson) match {
        case Success(JsObject(inputs)) => validateLabelRestrictions(inputs).map(_ => labels)
        case Failure(reason: Throwable) => s"Workflow contains invalid labels JSON: ${reason.getMessage}".invalidNel
        case _ => """Invalid workflow labels JSON. Expected a JsObject of "labelKey": "labelValue" values.""".invalidNel
      }
    }

    partialSources match {
      case Valid(partialSource) =>
        (validateInputs(partialSource),
          validateOptions(partialSource.workflowOptions),
          validateLabels(partialSource.customLabels.getOrElse("{}")),
          partialSource.workflowUrl.traverse(validateWorkflowUrl)) mapN {
          case (wfInputs, wfOptions, workflowLabels, wfUrl) =>
            wfInputs.map(inputsJson => WorkflowSourceFilesCollection(
              workflowSource = partialSource.workflowSource,
              workflowUrl = wfUrl,
              workflowRoot = partialSource.workflowRoot,
              workflowType = partialSource.workflowType,
              workflowTypeVersion = partialSource.workflowTypeVersion,
              inputsJson = inputsJson,
              workflowOptionsJson = wfOptions.asPrettyJson,
              labelsJson = workflowLabels,
              importsFile = partialSource.zippedImports,
              warnings = partialSource.warnings,
              workflowOnHold = partialSource.workflowOnHold))
        }
      case Invalid(err) => err.invalid
    }
  }

  def mergeMaps(allInputs: Seq[Option[String]]): JsObject = {
    val convertToMap = allInputs.map(x => toMap(x))
    JsObject(convertToMap reduce (_ ++ _))
  }

  def validateWorkflowUrl(workflowUrl: String): ErrorOr[WorkflowUrl] = {
    def convertStringToUrl(workflowUrl: String): ErrorOr[WorkflowUrl] = {
      Try(new URL(workflowUrl)) match {
        case Success(_) => workflowUrl.validNel
        case Failure(e) => s"Error while validating workflow url: ${e.getMessage}".invalidNel
      }
    }

    val len = workflowUrl.length
    if (len > MaxWorkflowUrlLength) s"Invalid workflow url: url has length $len, longer than the maximum allowed $MaxWorkflowUrlLength characters".invalidNel
    else convertStringToUrl(workflowUrl)
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

