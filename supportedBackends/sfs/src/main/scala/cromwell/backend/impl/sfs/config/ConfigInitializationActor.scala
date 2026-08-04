package cromwell.backend.impl.sfs.config

import cats.instances.future._
import cats.instances.option._
import cats.syntax.traverse._
import com.google.auth.oauth2.OAuth2Credentials
import common.validation.Validation._
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.sfs._
import cromwell.backend.standard.{
  StandardInitializationActorParams,
  StandardInitializationData,
  StandardValidatedRuntimeAttributesBuilder
}
import cromwell.cloudsupport.gcp.GoogleConfiguration
import net.ceedubs.ficus.Ficus._
import wdl.draft2.model.WdlNamespace

import scala.concurrent.Future

/**
  * Extension of the SharedFileSystemBackendInitializationData with declarations of extra runtime attributes, and a
  * wdl namespace containing various tasks for submitting, killing, etc.
  *
  * @param workflowPaths            The paths for the workflow.
  * @param runtimeAttributesBuilder The customized runtime attributes builder with extra validations for the
  *                                 declarations.
  * @param declarationValidations   A collection of validations for each declaration.
  * @param wdlNamespace             A collection of WDL tasks for submitting, killing, etc.
  */
class ConfigInitializationData(workflowPaths: WorkflowPaths,
                               runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder,
                               val googleRegistryCredentialsOption: Option[OAuth2Credentials],
                               val declarationValidations: Seq[DeclarationValidation],
                               val wdlNamespace: WdlNamespace
) extends StandardInitializationData(workflowPaths,
                                     runtimeAttributesBuilder,
                                     classOf[SharedFileSystemExpressionFunctions]
    )

/**
  * Extends the SharedFileSystemInitializationActor to create an instance of the ConfigInitializationData.
  *
  * The runtime attributes validation will be appended with validations based on the declarations from the config.
  *
  * @param params Parameters to create an initialization actor.
  */
class ConfigInitializationActor(params: StandardInitializationActorParams)
    extends SharedFileSystemInitializationActor(params) {

  private lazy val configWdlNamespace = new ConfigWdlNamespace(params.configurationDescriptor.backendConfig)

  lazy val declarationValidations: Seq[DeclarationValidation] =
    DeclarationValidation.fromDeclarations(configWdlNamespace.runtimeDeclarations,
                                           configWdlNamespace.callCachedRuntimeAttributes
    )

  /**
    * Return optional credentials for call caching GAR/GCR images.
    */
  private lazy val googleRegistryCredentialsOption: Future[Option[OAuth2Credentials]] = {
    val dockerGoogleAuthOption =
      standardParams.configurationDescriptor.backendConfig.getAs[String]("docker.google.auth")
    dockerGoogleAuthOption traverse { dockerGoogleAuth =>
      val googleConfiguration = GoogleConfiguration(standardParams.configurationDescriptor.globalConfig)
      val googleAuthTry =
        googleConfiguration
          .auth(dockerGoogleAuth)
          .toTry(s"Error retrieving google auth mode $dockerGoogleAuth")
      for {
        googleAuth <- Future.fromTry(googleAuthTry)
        credentials <- Future(
          googleAuth.credentials(
            workflowDescriptor.workflowOptions.get(_).get,
            List("https://www.googleapis.com/auth/cloud-platform")
          )
        )
      } yield credentials
    }
  }

  override lazy val initializationData: Future[ConfigInitializationData] = {
    val wdlNamespace = configWdlNamespace.wdlNamespace
    for {
      workflowPathsActual <- workflowPaths
      googleRegistryCredentialsOptionActual <- googleRegistryCredentialsOption
    } yield new ConfigInitializationData(
      workflowPathsActual,
      runtimeAttributesBuilder,
      googleRegistryCredentialsOptionActual,
      declarationValidations,
      wdlNamespace
    )
  }

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder = {
    val declared = declarationValidations.map(_.makeValidation())
    super.runtimeAttributesBuilder.withValidation(declared: _*)
  }
}
