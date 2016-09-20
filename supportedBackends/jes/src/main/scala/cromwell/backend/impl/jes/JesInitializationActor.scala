package cromwell.backend.impl.jes

import java.io.IOException

import akka.actor.{ActorRef, Props}
import com.google.api.services.genomics.Genomics
import cromwell.backend.impl.jes.JesInitializationActor._
import cromwell.backend.impl.jes.authentication.{GcsLocalizing, JesAuthInformation, JesCredentials}
import cromwell.backend.impl.jes.io._
import cromwell.backend.validation.RuntimeAttributesDefault
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.{BackendInitializationData, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.WorkflowOptions
import cromwell.core.retry.Retry
import cromwell.filesystems.gcs.{ClientSecrets, GoogleAuthMode}
import spray.json.JsObject
import wdl4s.types.{WdlBooleanType, WdlFloatType, WdlIntegerType, WdlStringType}
import wdl4s.values.WdlValue
import wdl4s.{Call, WdlExpression}

import scala.concurrent.Future
import scala.util.Try

object JesInitializationActor {
  val SupportedKeys = Set(CpuKey, MemoryKey, DockerKey, FailOnStderrKey, ContinueOnReturnCodeKey, JesRuntimeAttributes.ZonesKey,
    JesRuntimeAttributes.PreemptibleKey, JesRuntimeAttributes.BootDiskSizeKey, JesRuntimeAttributes.DisksKey)

  def props(workflowDescriptor: BackendWorkflowDescriptor,
            calls: Seq[Call],
            jesConfiguration: JesConfiguration,
            serviceRegistryActor: ActorRef): Props =
    Props(new JesInitializationActor(workflowDescriptor, calls, jesConfiguration, serviceRegistryActor: ActorRef))
}

class JesInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
                             override val calls: Seq[Call],
                             private[jes] val jesConfiguration: JesConfiguration,
                             override val serviceRegistryActor: ActorRef)
  extends BackendWorkflowInitializationActor {

  override protected def runtimeAttributeValidators: Map[String, (Option[WdlExpression]) => Boolean] = Map(
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
      clientSecrets <- List(jesConfiguration.jesAttributes.gcsFilesystemAuth) collectFirst { case s: ClientSecrets => s }
      token <- workflowDescriptor.workflowOptions.get(GoogleAuthMode.RefreshTokenOptionKey).toOption
    } yield GcsLocalizing(clientSecrets, token)
  }

  private val iOExecutionContext = context.system.dispatchers.lookup(IoDispatcher)


  override protected def coerceDefaultRuntimeAttributes(options: WorkflowOptions): Try[Map[String, WdlValue]] = {
    RuntimeAttributesDefault.workflowOptionsDefault(options, JesRuntimeAttributes.coercionMap)
  }

  /**
    * A call which happens before anything else runs
    */
  override def beforeAll(): Future[Option[BackendInitializationData]] = {

    val genomicsCredential = jesConfiguration.jesAttributes.genomicsCredential(workflowDescriptor.workflowOptions)
    val gcsCredential = jesConfiguration.jesAttributes.gcsCredential(workflowDescriptor.workflowOptions)

    val jesCredentials = JesCredentials(genomicsCredential = genomicsCredential, gcsCredential = gcsCredential)
    def buildGenomics: Future[Genomics] = Future {
      GenomicsFactory(jesConfiguration.googleConfig.applicationName, genomicsCredential, jesConfiguration.jesAttributes.endpointUrl)
    }

    for {
      // generate single filesystem and genomics instances
      genomics <- buildGenomics
      workflowPaths = new JesWorkflowPaths(workflowDescriptor, jesConfiguration, jesCredentials)(iOExecutionContext)
      _ <- if (jesConfiguration.needAuthFileUpload) writeAuthenticationFile(workflowPaths) else Future.successful(())
      _ = publishWorkflowRoot(workflowPaths.workflowRootPath.toString)
    } yield Option(JesBackendInitializationData(workflowPaths, genomics))
  }

  private def writeAuthenticationFile(workflowPath: JesWorkflowPaths): Future[Unit] = {
    generateAuthJson(jesConfiguration.dockerCredentials, refreshTokenAuth) map { content =>
      val path = workflowPath.gcsAuthFilePath
      val upload = () => Future(path.writeAsJson(content))

      workflowLogger.info(s"Creating authentication file for workflow ${workflowDescriptor.id} at \n ${path.toString}")
      Retry.withRetry(upload, isFatal = isFatalJesException, isTransient = isTransientJesException)(context.system) map { _ => () } recover {
        case failure => throw new IOException("Failed to upload authentication file", failure)
      }
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
