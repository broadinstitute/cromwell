package cromwell.backend.impl.jes

import akka.actor.Props
import com.google.api.services.genomics.Genomics
import cromwell.backend.impl.jes.JesImplicits.GoogleAuthWorkflowOptions
import cromwell.backend.impl.jes.JesInitializationActor._
import cromwell.backend.impl.jes.authentication.{GcsLocalizing, JesAuthInformation}
import cromwell.backend.impl.jes.io._
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.{BackendInitializationData, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.core.retry.Retry
import cromwell.filesystems.gcs.{ClientSecrets, GcsFileSystem, GcsFileSystemProvider, GoogleAuthMode}
import spray.json.JsObject
import wdl4s.types.{WdlBooleanType, WdlFloatType, WdlIntegerType, WdlStringType}
import wdl4s.{Call, WdlExpression}

import scala.concurrent.Future

object JesInitializationActor {
  val SupportedKeys = Set(CpuKey, MemoryKey, DockerKey, FailOnStderrKey, ContinueOnReturnCodeKey, JesRuntimeAttributes.ZonesKey,
    JesRuntimeAttributes.PreemptibleKey, JesRuntimeAttributes.BootDiskSizeKey, JesRuntimeAttributes.DisksKey)

  def props(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], jesConfiguration: JesConfiguration): Props =
    Props(new JesInitializationActor(workflowDescriptor, calls, jesConfiguration)).withDispatcher("akka.dispatchers.slow-actor-dispatcher")
}

class JesInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
                             override val calls: Seq[Call],
                             private[jes] val jesConfiguration: JesConfiguration)
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

  //TODO PBE: Workflow options may need to be validated for JES.

  /**
    * A call which happens before anything else runs
    */
  override def beforeAll(): Future[Option[BackendInitializationData]] = {

    def buildGenomics: Future[Genomics] = Future {
      val jesAttributes = jesConfiguration.jesAttributes
      GenomicsFactory(
        jesConfiguration.googleConfig, jesAttributes.genomicsAuth.credential(workflowDescriptor.workflowOptions.toGoogleAuthOptions), jesAttributes.endpointUrl)
    }

    def buildGcsFileSystem: Future[GcsFileSystem] = Future {
      val storage = jesConfiguration.jesAttributes.gcsFilesystemAuth.buildStorage(
        workflowDescriptor.workflowOptions.toGoogleAuthOptions, jesConfiguration.googleConfig)
      GcsFileSystem(GcsFileSystemProvider(storage))
    }

    for {
      // generate single filesystem and genomics instances
      gcsFileSystem <- buildGcsFileSystem
      genomics <- buildGenomics
      workflowPaths = new JesWorkflowPaths(workflowDescriptor, jesConfiguration, gcsFileSystem)
      _ <- if (jesConfiguration.needAuthFileUpload) writeAuthenticationFile(workflowPaths) else Future.successful(())
      _ = publishWorkflowRoot(workflowPaths.workflowRootPath.toString)
    } yield Option(JesBackendInitializationData(workflowPaths, genomics))
  }

  private def writeAuthenticationFile(workflowPath: JesWorkflowPaths): Future[Unit] = {
    generateAuthJson(jesConfiguration.dockerCredentials, refreshTokenAuth) map { content =>
      val path = workflowPath.gcsAuthFilePath
      val upload = () => Future(path.writeAsJson(content))

      workflowLogger.info(s"Creating authentication file for workflow ${workflowDescriptor.id} at \n ${path.toString}")
      Retry.withRetry(upload, isFatal = isFatalJesException, isTransient = isTransientJesException)(context.system) map { _ => () }
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
