package cromwell.backend.impl.jes

import java.io.IOException

import akka.actor.{ActorRef, Props}
import cats.instances.future._
import cats.syntax.functor._
import com.google.api.services.genomics.Genomics
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.backend.impl.jes.JesInitializationActor._
import cromwell.backend.impl.jes.authentication.{GcsLocalizing, JesAuthInformation}
import cromwell.backend.impl.jes.io._
import cromwell.backend.validation.RuntimeAttributesDefault
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.{BackendInitializationData, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.core.WorkflowOptions
import cromwell.filesystems.gcs.auth.{ClientSecrets, GoogleAuthMode}
import spray.json.JsObject
import wdl4s.TaskCall
import wdl4s.types.{WdlBooleanType, WdlFloatType, WdlIntegerType, WdlStringType}
import wdl4s.values.WdlValue

import scala.concurrent.Future
import scala.language.postfixOps
import scala.util.Try

object JesInitializationActor {
  val SupportedKeys = Set(CpuKey, MemoryKey, DockerKey, FailOnStderrKey, ContinueOnReturnCodeKey, JesRuntimeAttributes.ZonesKey,
    JesRuntimeAttributes.PreemptibleKey, JesRuntimeAttributes.BootDiskSizeKey, JesRuntimeAttributes.DisksKey)

  def props(workflowDescriptor: BackendWorkflowDescriptor,
            calls: Set[TaskCall],
            jesConfiguration: JesConfiguration,
            serviceRegistryActor: ActorRef): Props =
    Props(new JesInitializationActor(workflowDescriptor, calls, jesConfiguration, serviceRegistryActor: ActorRef)).withDispatcher(BackendDispatcher)
}

class JesInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
                             override val calls: Set[TaskCall],
                             private[jes] val jesConfiguration: JesConfiguration,
                             override val serviceRegistryActor: ActorRef)
  extends BackendWorkflowInitializationActor {

  override protected def runtimeAttributeValidators: Map[String, (Option[WdlValue]) => Boolean] = Map(
    CpuKey -> wdlTypePredicate(valueRequired = false, WdlIntegerType.isCoerceableFrom),
    MemoryKey -> wdlTypePredicate(valueRequired = false, WdlStringType.isCoerceableFrom),
    DockerKey -> wdlTypePredicate(valueRequired = true, WdlStringType.isCoerceableFrom),
    FailOnStderrKey -> wdlTypePredicate(valueRequired = false, WdlBooleanType.isCoerceableFrom),
    ContinueOnReturnCodeKey -> continueOnReturnCodePredicate(valueRequired = false),
    JesRuntimeAttributes.PreemptibleKey -> wdlTypePredicate(valueRequired = false, WdlIntegerType.isCoerceableFrom),
    JesRuntimeAttributes.BootDiskSizeKey -> wdlTypePredicate(valueRequired = false, WdlFloatType.isCoerceableFrom),

    // TODO (eventually): make these more appropriate pre-checks
    JesRuntimeAttributes.ZonesKey -> wdlTypePredicate(valueRequired = false, WdlStringType.isCoerceableFrom),
    JesRuntimeAttributes.DisksKey -> wdlTypePredicate(valueRequired = false, WdlStringType.isCoerceableFrom))

  override val configurationDescriptor = jesConfiguration.configurationDescriptor

  private[jes] lazy val refreshTokenAuth: Option[JesAuthInformation] = {
    for {
      clientSecrets <- List(jesConfiguration.jesAttributes.auths.gcs) collectFirst { case s: ClientSecrets => s }
      token <- workflowDescriptor.workflowOptions.get(GoogleAuthMode.RefreshTokenOptionKey).toOption
    } yield GcsLocalizing(clientSecrets, token)
  }

  override protected def coerceDefaultRuntimeAttributes(options: WorkflowOptions): Try[Map[String, WdlValue]] = {
    RuntimeAttributesDefault.workflowOptionsDefault(options, JesRuntimeAttributes.coercionMap)
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
    } yield Option(JesBackendInitializationData(workflowPaths, genomics))
  }

  private def writeAuthenticationFile(workflowPath: JesWorkflowPaths): Future[Unit] = {
    generateAuthJson(jesConfiguration.dockerCredentials, refreshTokenAuth) map { content =>
      val path = workflowPath.gcsAuthFilePath
      workflowLogger.info(s"Creating authentication file for workflow ${workflowDescriptor.id} at \n ${path.toUri}")
      Future(path.writeAsJson(content)).void.recoverWith {
        case failure => Future.failed(new IOException("Failed to upload authentication file", failure))
      } void
    } getOrElse Future.successful(())
  }

  def generateAuthJson(authInformation: Option[JesAuthInformation]*) = {
    authInformation.flatten map { _.toMap } match {
      case Nil => None
      case jsons =>
        val authsValues = jsons.reduce(_ ++ _) mapValues JsObject.apply
        Option(JsObject("auths" -> JsObject(authsValues)).prettyPrint)
    }
  }

  /**
    * Validate that this WorkflowBackendActor can run all of the calls that it's been assigned
    */
  override def validate(): Future[Unit] = {
    Future {
      calls foreach { call =>
        val runtimeAttributes = call.task.runtimeAttributes.attrs
        val notSupportedAttributes = runtimeAttributes filterKeys { !SupportedKeys.contains(_) }

        if (notSupportedAttributes.nonEmpty) {
          val notSupportedAttrString = notSupportedAttributes.keys mkString ", "
          workflowLogger.warn(s"Key/s [$notSupportedAttrString] is/are not supported by JesBackend. Unsupported attributes will not be part of jobs executions.")
        }
      }
    }
  }
}
