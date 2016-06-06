package cromwell.backend.impl.jes

import akka.actor.Props
import cromwell.backend.impl.jes.JesInitializationActor._
import cromwell.backend.impl.jes.authentication.{GcsLocalizing, JesAuthInformation}
import cromwell.backend.impl.jes.io._
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.{BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.core.retry.Retry
import cromwell.filesystems.gcs.{ClientSecrets, GoogleAuthMode}
import spray.json.JsObject
import wdl4s.Call

import scala.concurrent.Future

object JesInitializationActor {
  val SupportedKeys = Set(CpuKey, MemoryKey, DockerKey, FailOnStderrKey, ContinueOnReturnCodeKey, JesRuntimeAttributes.ZonesKey,
    JesRuntimeAttributes.PreemptibleKey, JesRuntimeAttributes.BootDiskSizeKey, JesRuntimeAttributes.DisksKey)

  def props(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], jesConfiguration: JesConfiguration): Props =
    Props(new JesInitializationActor(workflowDescriptor, calls, jesConfiguration))
}

class JesInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
                             override val calls: Seq[Call],
                             jesConfiguration: JesConfiguration) extends BackendWorkflowInitializationActor {

  override val configurationDescriptor = jesConfiguration.configurationDescriptor
  private val workflowPaths = new JesWorkflowPaths(workflowDescriptor, jesConfiguration)

  private lazy val refreshTokenAuth: Option[JesAuthInformation] = {
    for {
      clientSecrets <- List(jesConfiguration.jesAttributes.gcsFilesystemAuth) collectFirst { case s: ClientSecrets => s }
      token <- workflowDescriptor.workflowOptions.get(GoogleAuthMode.RefreshTokenOptionKey).toOption
    } yield GcsLocalizing(clientSecrets, token)
  }

  /**
    * Abort all initializations.
    */
  override def abortInitialization(): Unit = ???

  //TODO: Workflow options may need to be validated for JES.

  /**
    * A call which happens before anything else runs
    */
  override def beforeAll(): Future[Unit] = {
    if (jesConfiguration.needAuthFileUpload) writeAuthenticationFile()
    else Future.successful(())
  }

  private def writeAuthenticationFile(): Future[Unit] = {
    generateAuthJson(jesConfiguration.dockerCredentials, refreshTokenAuth) map { content =>
      val path = workflowPaths.gcsAuthFilePath
      val upload = () => Future(path.writeAsJson(content))

      log.info(s"Creating authentication file for workflow ${workflowDescriptor.id} at \n ${path.toString}")
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
          log.warning(s"Key/s [$notSupportedAttrString] is/are not supported by JesBackend. Unsupported attributes will not be part of jobs executions.")
        }

        runtimeAttributes.get(DockerKey).orElse(throw new IllegalArgumentException(s"$DockerKey mandatory runtime attribute is missing."))
      }
    }
  }
}
