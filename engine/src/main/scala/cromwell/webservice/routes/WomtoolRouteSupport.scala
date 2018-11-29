package cromwell.webservice.routes

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model._
import akka.http.scaladsl.server.Directives._
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import cromwell.core.WorkflowSourceFilesCollection
import cromwell.engine.language.CromwellLanguages
import cromwell.languages.LanguageFactory
import cromwell.webservice.WebServiceUtils
import cromwell.webservice.WebServiceUtils.EnhancedThrowable
import spray.json.{JsArray, JsBoolean, JsObject, JsString}
import wom.core.WorkflowSource

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

trait WomtoolRouteSupport extends WebServiceUtils {

  implicit val ec: ExecutionContext
  implicit val materializer: ActorMaterializer
  implicit val timeout: Timeout

  val womtoolRoutes =
    path("womtool" / Segment / "describe") { _ =>
      post {
        entity(as[Multipart.FormData]) { formData: Multipart.FormData =>
          onComplete(materializeFormData(formData)) {
            case Success(data) =>

              // TODO: move constants to WebServiceUtils, adopt in PartialWorkflowSources
              val workflowSource = data.get("workflowSource").map(_.utf8String)
              val workflowUrl = data.get("workflowUrl").map(_.utf8String)
              val workflowInputs = data.get("workflowInputs").map(_.utf8String)
              val workflowType = data.get("workflowType").map(_.utf8String)
              val workflowVersion = data.get("workflowTypeVersion").map(_.utf8String)

              val wsfc = WorkflowSourceFilesCollection(
                workflowSource,
                workflowUrl,
                workflowRoot = None,
                workflowType,
                workflowVersion,
                workflowInputs.getOrElse(""),
                workflowOptionsJson = "",
                labelsJson = "",
                importsFile = None,
                workflowOnHold = false,
                warnings = Seq.empty
              )

              // TODO: would be nice to reuse findWorkflowSource
              (workflowSource, workflowUrl) match {
                case (Some(_), Some(_)) =>
                  new IllegalArgumentException("Must submit exactly one of workflow source, workflow URL; received both").failRequest(StatusCodes.BadRequest)
                case (None, Some(_)) =>
                  new Exception("URL submissions not yet supported").failRequest(StatusCodes.NotImplemented)
                case (Some(source), None) =>
                  complete {
                    // After much discussion, it was resolved that returning Future(something) is the right way to do
                    // things in Akka HTTP, when no legacy or migration constraints funnel one into using PerRequest
                    doExpensiveWDLValidation(source, wsfc) map { case (valid: Boolean, errors: List[String]) =>
                      JsObject(
                        Map(
                          "valid" -> JsBoolean(valid),
                          "errors" -> JsArray(errors.map(JsString(_)).toVector)
                        )
                      )
                    }
                  }
                case (None, None) =>
                  new IllegalArgumentException("Must submit exactly one of workflow source, workflow URL; received neither").failRequest(StatusCodes.BadRequest)
              }
            case Failure(e) => e.failRequest(StatusCodes.InternalServerError)
          }
        }
      }
    }

  private def doExpensiveWDLValidation(workflow: WorkflowSource, workflowSourceFilesCollection: WorkflowSourceFilesCollection): Future[(Boolean, List[String])] = {
    chooseFactory(workflow, workflowSourceFilesCollection) match {
      case Valid(factory: LanguageFactory) =>
        // Why do we pass in the rest of the language factories here? I cannot figure out what we ever use them for.
        Future(factory.getWomBundle(workflow, "{}", List.empty, List.empty).isRight)

        Future {
          factory.getWomBundle(workflow, "{}", List.empty, List.empty) match {
            case Right(_) => (true, List.empty)
            case Left(errors) => (false, errors.toList)
          }
        }
      case Invalid(e) =>
        Future.failed(new Exception(e.toList.mkString(", ")))
    }
  }

  // TODO: refactor me into somewhere central! (see also MaterializeWorkflowDescriptorActor)
  private def chooseFactory(workflowSource: WorkflowSource, wsfc: WorkflowSourceFilesCollection): ErrorOr[LanguageFactory] = {
    wsfc.workflowType match {
      case Some(languageName) if CromwellLanguages.instance.languages.contains(languageName.toUpperCase) =>
        val language = CromwellLanguages.instance.languages(languageName.toUpperCase)
        wsfc.workflowTypeVersion match {
          case Some(v) if language.allVersions.contains(v) => language.allVersions(v).valid
          case Some(other) => s"Unknown version '$other' for workflow language '$languageName'".invalidNel
          case _ => chooseFactory(workflowSource, wsfc).getOrElse(language.default).valid
        }
      case Some(other) => s"Unknown workflow type: $other".invalidNel[LanguageFactory]
      case None =>
        val allFactories = CromwellLanguages.instance.languages.values.flatMap(_.allVersions.values)
        allFactories.find(_.looksParsable(workflowSource)).getOrElse(CromwellLanguages.instance.default.default).validNel
    }
  }

}
