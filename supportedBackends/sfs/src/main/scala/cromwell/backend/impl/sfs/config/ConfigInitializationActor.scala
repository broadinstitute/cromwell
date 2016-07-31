package cromwell.backend.impl.sfs.config

import cromwell.backend.io.WorkflowPaths
import cromwell.backend.sfs._
import wdl4s.WdlNamespace

/**
  * Extension of the SharedFileSystemBackendInitializationData with declarations of extra runtime attributes, and a
  * wdl namespace containing various tasks for submiting, killing, etc.
  *
  * @param workflowPaths            The paths for the workflow.
  * @param runtimeAttributesBuilder The customized runtime attributes builder with extra validations for the
  *                                 declarations.
  * @param declarationValidations   A collection of validations for each declaration.
  * @param wdlNamespace             A collection of WDL tasks for submitting, killing, etc.
  */
class ConfigInitializationData
(
  workflowPaths: WorkflowPaths,
  runtimeAttributesBuilder: SharedFileSystemValidatedRuntimeAttributesBuilder,
  val declarationValidations: Seq[DeclarationValidation],
  val wdlNamespace: WdlNamespace)
  extends SharedFileSystemBackendInitializationData(workflowPaths, runtimeAttributesBuilder)

/**
  * Extends the SharedFileSystemInitializationActor to create an instance of the ConfigInitializationData.
  *
  * The runtime attributes validation will be appended with validations based on the declarations from the config.
  *
  * @param params Parameters to create an initialization actor.
  */
class ConfigInitializationActor(params: SharedFileSystemInitializationActorParams)
  extends SharedFileSystemInitializationActor(params) {

  lazy val configWdlNamespace = new ConfigWdlNamespace(params.configurationDescriptor.backendConfig)

  lazy val declarationValidations: Seq[DeclarationValidation] = {
    DeclarationValidation.fromDeclarations(configWdlNamespace.runtimeDeclarations)
  }

  override lazy val initializationData = {
    val wdlNamespace = configWdlNamespace.wdlNamespace
    new ConfigInitializationData(workflowPaths, runtimeAttributesBuilder, declarationValidations, wdlNamespace)
  }

  override lazy val runtimeAttributesBuilder = {
    val declared = declarationValidations.map(_.makeValidation())
    super.runtimeAttributesBuilder.withValidation(declared: _*)
  }
}
