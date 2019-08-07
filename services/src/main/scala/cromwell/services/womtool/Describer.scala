package cromwell.services.womtool

import cats.data.Validated.{Invalid, Valid}
import cromwell.core.WorkflowSourceFilesCollection
import cromwell.languages.util.ImportResolver.HttpResolver
import cromwell.languages.util.{ImportResolver, LanguageFactoryUtil}
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import cromwell.services.womtool.WomtoolServiceMessages.{DescribeFailure, DescribeResult, DescribeSuccess}
import cromwell.services.womtool.models.WorkflowDescription
import wom.core.WorkflowSource
import wom.executable.WomBundle
import wom.expression.NoIoFunctionSet

object Describer {

  def describeWorkflow(wsfc: WorkflowSourceFilesCollection): DescribeResult = {

    // The HTTP resolver is used to pull down workflows submitted by URL
    LanguageFactoryUtil.findWorkflowSource(wsfc.workflowSource, wsfc.workflowUrl, List(HttpResolver(None, Map.empty))) match {
      case Right(sourceAndResolvers: (WorkflowSource, List[ImportResolver.ImportResolver])) =>
        LanguageFactoryUtil.chooseFactory(sourceAndResolvers._1, wsfc) match {
          case Valid(factory: LanguageFactory) =>
            DescribeSuccess(
              description = describeWorkflowInner(factory, sourceAndResolvers._1, wsfc)
            )
          case Invalid(e) =>
            // `chooseFactory` generally should not fail, because we still choose the default factory even if `looksParsable`
            // returns uniformly `false`. (If the user submits gibberish, we will use the default factory and fail elsewhere.)
            // We could get here if there's a problem loading the factories' classes or the language configuration, however.
            // It is a BadRequest because ultimately the user submitted something we could not understand.
            DescribeFailure(
              reason = e.toList.mkString(", ")
            )
        }
      case Left(errors) =>
        // Likely reasons: no workflow provided, workflow URL could not be read
        DescribeFailure(
          reason = errors.toList.mkString(", ")
        )
    }
  }

  // By this point there are no "out of band" errors that can occur (i.e. those that would indicate a BadRequest, versus just showing up in the `errors` list)
  private def describeWorkflowInner(factory: LanguageFactory, workflow: WorkflowSource, workflowSourceFilesCollection: WorkflowSourceFilesCollection): WorkflowDescription = {

    // Mirror of the inputs/no inputs fork in womtool.validate.Validate
    if (workflowSourceFilesCollection.inputsJson.isEmpty) {
      // No inputs: just load up the WomBundle
      factory.getWomBundle(workflow, workflowSourceOrigin = None, workflowOptionsJson = "{}", List(HttpResolver(None, Map.empty)), List(factory)) match {
        case Right(bundle: WomBundle) =>
          WorkflowDescription.fromBundle(bundle, factory.languageName, factory.languageVersionName, List.empty)
        case Left(workflowErrors) =>
          WorkflowDescription.withErrors(workflowErrors.toList)
      }
    } else {
      // Inputs: load up the WomBundle and then try creating an executable with WomBundle + inputs
      factory.getWomBundle(workflow, workflowSourceOrigin = None, workflowOptionsJson = "{}", List(HttpResolver(None, Map.empty)), List(factory)) match {
        case Right(bundle) =>
          factory.createExecutable(bundle, workflowSourceFilesCollection.inputsJson, NoIoFunctionSet) match {
            // Throw away the executable, all we care about is whether it created successfully (i.e. the inputs are valid)
            case Right(_: ValidatedWomNamespace) =>
              WorkflowDescription.fromBundle(bundle, factory.languageName, factory.languageVersionName)
            case Left(inputErrors) =>
              WorkflowDescription.fromBundle(bundle, factory.languageName, factory.languageVersionName, inputErrors.toList)
          }
        case Left(workflowErrors) =>
          WorkflowDescription.withErrors(workflowErrors.toList)
      }
    }

  }

}
