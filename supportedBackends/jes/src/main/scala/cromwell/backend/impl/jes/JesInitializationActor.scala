package cromwell.backend.impl.jes

import java.io.IOException

import akka.actor.{ActorRef, Props}
import cats.instances.future._
import cats.syntax.functor._
import com.google.api.services.genomics.Genomics
import cromwell.backend.impl.jes.authentication.{GcsLocalizing, JesAuthInformation}
import cromwell.backend.impl.jes.io._
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.filesystems.gcs.auth.{ClientSecrets, GoogleAuthMode}
import spray.json.JsObject
import wdl4s.TaskCall

import scala.concurrent.Future
import scala.language.postfixOps

object JesInitializationActor {
  /* NOTE: Only used by tests */
  def props(workflowDescriptor: BackendWorkflowDescriptor,
            calls: Set[TaskCall],
            jesConfiguration: JesConfiguration,
            serviceRegistryActor: ActorRef): Props = {
    val params = JesInitializationActorParams(workflowDescriptor, calls, jesConfiguration, serviceRegistryActor)
    Props(new JesInitializationActor(params)).withDispatcher(BackendDispatcher)
  }
}

case class JesInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  calls: Set[TaskCall],
  jesConfiguration: JesConfiguration,
  serviceRegistryActor: ActorRef
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = jesConfiguration.configurationDescriptor
}

class JesInitializationActor(jesParams: JesInitializationActorParams)
  extends StandardInitializationActor with JesAsyncIo {

  override val standardParams: StandardInitializationActorParams = jesParams

  private val jesConfiguration = jesParams.jesConfiguration

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    JesRuntimeAttributes.runtimeAttributesBuilder(jesConfiguration)

  private[jes] lazy val refreshTokenAuth: Option[JesAuthInformation] = {
    for {
      clientSecrets <- List(jesConfiguration.jesAttributes.auths.gcs) collectFirst { case s: ClientSecrets => s }
      token <- workflowDescriptor.workflowOptions.get(GoogleAuthMode.RefreshTokenOptionKey).toOption
    } yield GcsLocalizing(clientSecrets, token)
  }

  /**
    * A call which happens before anything else runs
    */
  override def beforeAll(): Future[Option[BackendInitializationData]] = {

    def buildGenomics: Future[Genomics] = Future {
      jesConfiguration.genomicsFactory.withOptions(workflowDescriptor.workflowOptions)
    }

    for {
      genomics <- buildGenomics
      workflowPaths = new JesWorkflowPaths(workflowDescriptor, jesConfiguration)(context.system)
      _ <- if (jesConfiguration.needAuthFileUpload) writeAuthenticationFile(workflowPaths) else Future.successful(())
      _ = publishWorkflowRoot(workflowPaths.workflowRoot.toString)
    } yield Option(JesBackendInitializationData(workflowPaths, runtimeAttributesBuilder, jesConfiguration, genomics))
  }

  private def writeAuthenticationFile(workflowPath: JesWorkflowPaths): Future[Unit] = {
    generateAuthJson(jesConfiguration.dockerCredentials, refreshTokenAuth) map { content =>
      val path = workflowPath.gcsAuthFilePath
      workflowLogger.info(s"Creating authentication file for workflow ${workflowDescriptor.id} at \n ${path.toUri}")
      writeAsJson(path, content).void.recoverWith {
        case failure => Future.failed(new IOException("Failed to upload authentication file", failure))
      } void
    } getOrElse Future.successful(())
  }

  def generateAuthJson(authInformation: Option[JesAuthInformation]*): Option[String] = {
    authInformation.flatten map { _.toMap } match {
      case Nil => None
      case jsons =>
        val authsValues = jsons.reduce(_ ++ _) mapValues JsObject.apply
        Option(JsObject("auths" -> JsObject(authsValues)).prettyPrint)
    }
  }
}
