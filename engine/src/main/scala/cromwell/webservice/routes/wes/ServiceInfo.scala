package cromwell.webservice.routes.wes

import java.net.URL

import akka.actor.ActorRef
import akka.util.Timeout
import cromwell.core.filesystem.CromwellFileSystems
import cromwell.languages.config.CromwellLanguages
import cromwell.webservice.routes.CromwellApiService
import net.ceedubs.ficus.Ficus._
import com.typesafe.config.ConfigFactory
import spray.json.DefaultJsonProtocol
import akka.pattern.ask
import cromwell.core.WorkflowState
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.GetWorkflowStoreStats
import cromwell.webservice.routes.wes.WesState.WesState

import scala.concurrent.{ExecutionContext, Future}

object ServiceInfo {
  private lazy val config = ConfigFactory.load()

  // Which languages and which versions of those languages does this Cromwell support?
  // CromwellLanguages.instance.languages is a Map[Language, Map[Version, Factory]]. We want Map[Language -> Version]
  lazy val WorkflowTypeVersion = CromwellLanguages.instance.languages.map(x => x._1 -> x._2.allVersions.keys)

  // Which versions of WES does this Cromwell support?
  val SupportedWesVersions = List("1.0")

  // Which filesystems does this Cromwell support? NB: It's a known flaw of spec that not all backends may support all FS
  lazy val SupportedFilesystemProtocols = CromwellFileSystems.instance.supportedFileSystems

  // In Cromwell terms, default workflow options. Viewing this as too much trouble to fill in unless someone complains
  val DefaultWorkflowEngineParameters = List.empty[DefaultWorkflowEngineParameter]

  // Which engines and which version(s) of those engines does this WES support. Kinda obvious for Cromwell
  val WorkflowEngineVerisons = Map("Cromwell" -> CromwellApiService.cromwellVersion)

  // URL which provides site specific authorization information. Doesn't really apply to Cromwell for now
  lazy val AuthInstructionsUrl = config.as[URL]("ga4gh.wes.auth-instructions-url")
  // URL which provides site specific contact information for security issues
  lazy val ContactInfoUrl = config.as[URL]("ga4gh.wes.contact-info-url")

  /*
    Optional key/value pairs. Could potentially be added to the config but not going to bother until if/when someone asks for it
   */
  val Tags = Map.empty[String, String]

  /**
    * Generate any runtime level information and create a response to the client
    */
  def toWesResponse(workflowStoreActor: ActorRef)(implicit ec: ExecutionContext, timeout: Timeout): Future[WesStatusInfoResponse] = {
    workflowStats(workflowStoreActor).map(stats =>
      WesStatusInfoResponse(WorkflowTypeVersion,
        SupportedWesVersions,
        SupportedFilesystemProtocols,
        WorkflowEngineVerisons,
        DefaultWorkflowEngineParameters,
        stats,
        AuthInstructionsUrl.toString,
        ContactInfoUrl.toString,
        Tags)
    )
  }

  /**
    * Retrieve a map from state to count for all represented non-terminal workflow states
    */
  private def workflowStats(workflowStoreActor: ActorRef)(implicit ec: ExecutionContext, timeout: Timeout): Future[Map[WesState, Int]] = {
    workflowStoreActor.ask(GetWorkflowStoreStats)
      .mapTo[Map[WorkflowState, Int]]
      .map(m => m.map(e => WesState.fromCromwellStatus(e._1) -> e._2)) // Convert WorkflowState -> WesState
  }
}

/*
  Models the type DefaultWorkflowEngineParameter in the WES Swagger
 */
final case class DefaultWorkflowEngineParameter(name: String, `type`: String, default_value: String)

object DefaultWorkflowEngineParameter extends DefaultJsonProtocol {
  implicit val DefaultWorkflowEngineParameterFormat = jsonFormat3(DefaultWorkflowEngineParameter.apply)
}
