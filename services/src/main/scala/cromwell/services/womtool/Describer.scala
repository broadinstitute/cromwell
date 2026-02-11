package cromwell.services.womtool

import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import cromwell.core.WorkflowSourceFilesCollection
import cromwell.languages.util.ImportResolver.{
  zippedImportResolver,
  DirectoryResolver,
  HttpResolver,
  ImportAuthProvider
}
import common.validation.ErrorOr.ErrorOr
import cromwell.languages.util.{ImportResolver, LanguageFactoryUtil}
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import cromwell.services.womtool.WomtoolServiceMessages.{DescribeFailure, DescribeResult, DescribeSuccess}
import cromwell.services.womtool.models.WorkflowDescription
import wom.core.WorkflowSource
import wom.executable.WomBundle
import wom.expression.NoIoFunctionSet

object Describer {

  def describeWorkflow(wsfc: WorkflowSourceFilesCollection, authProviders: List[ImportAuthProvider]): DescribeResult = {

    val zippedImportResolverEither: ErrorOr[Option[DirectoryResolver]] =
      wsfc.importsZipFileOption match {
        case None => None.validNel
        case Some(zipContent) => zippedImportResolver(zipContent).map(Option.apply)
      }

    zippedImportResolverEither match {
      case Valid(zippedImportResolverOpt) =>
        val importResolvers = zippedImportResolverOpt.toList :+ HttpResolver(
          None,
          Map.empty,
          authProviders
        )
        LanguageFactoryUtil.findWorkflowSource(wsfc.workflowSource, wsfc.workflowUrl, importResolvers) match {
          case Right((workflowSource: WorkflowSource, importResolvers: List[ImportResolver.ImportResolver])) =>
            LanguageFactoryUtil.chooseFactory(workflowSource, wsfc) match {
              case Valid(factory: LanguageFactory) =>
                DescribeSuccess(
                  description = describeWorkflowInner(factory, workflowSource, importResolvers, wsfc)
                )
              // `chooseFactory` generally should not fail, because we still choose the default factory even if `looksParsable`
              // returns uniformly `false`. (If the user submits gibberish, we will use the default factory and fail elsewhere.)
              // We could get here if there's a problem loading the factories' classes or the language configuration, however.
              // It is a BadRequest because ultimately the user submitted something we could not understand.
              case Invalid(e) =>
                DescribeFailure(
                  reason = e.toList.mkString(", ")
                )
            }
          // Likely reasons: no workflow provided, workflow URL could not be read
          case Left(errors) =>
            DescribeFailure(
              reason = errors.toList.mkString(", ")
            )
        }
      case Invalid(errors) =>
        DescribeFailure(
          reason = errors.toList.mkString(", ")
        )
    }
  }

  // By this point there are no "out of band" errors that can occur (i.e. those that would indicate a BadRequest, versus just showing up in the `errors` list)
  private def describeWorkflowInner(factory: LanguageFactory,
                                    workflowSource: WorkflowSource,
                                    importResolvers: List[ImportResolver.ImportResolver],
                                    workflowSourceFilesCollection: WorkflowSourceFilesCollection
  ): WorkflowDescription = {

    val submittedDescriptorType = Map(
      "descriptorType" -> factory.languageName,
      "descriptorTypeVersion" -> factory.languageVersionName
    )

    // Mirror of the inputs/no inputs fork in womtool.validate.Validate
    if (workflowSourceFilesCollection.inputsJson.isEmpty) {
      // No inputs: just load up the WomBundle
      factory.getWomBundle(workflowSource,
                           workflowSourceOrigin = None,
                           workflowOptionsJson = "{}",
                           importResolvers,
                           List(factory)
      ) match {
        case Right(bundle: WomBundle) =>
          WorkflowDescription.fromBundle(bundle, submittedDescriptorType, List.empty)
        case Left(workflowErrors) =>
          WorkflowDescription.withErrors(workflowErrors.toList, submittedDescriptorType)
      }
    } else {
      // Inputs: load up the WomBundle and then try creating an executable with WomBundle + inputs
      factory.getWomBundle(workflowSource,
                           workflowSourceOrigin = None,
                           workflowOptionsJson = "{}",
                           importResolvers,
                           List(factory)
      ) match {
        case Right(bundle) =>
          factory.createExecutable(bundle, workflowSourceFilesCollection.inputsJson, NoIoFunctionSet) match {
            // Throw away the executable, all we care about is whether it created successfully (i.e. the inputs are valid)
            case Right(_: ValidatedWomNamespace) =>
              WorkflowDescription.fromBundle(bundle, submittedDescriptorType)
            case Left(inputErrors) =>
              WorkflowDescription.fromBundle(bundle, submittedDescriptorType, inputErrors.toList)
          }
        case Left(workflowErrors) =>
          WorkflowDescription.withErrors(workflowErrors.toList, submittedDescriptorType)
      }
    }

  }

}
