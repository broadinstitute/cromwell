package languages.wdl.draft2

// TODO Scala 2.11: cats.syntax.either._ again
import cats.syntax.either._
import cats.instances.either._
import cats.instances.list._
import cats.syntax.functor._
import cats.data.EitherT.fromEither
import cats.effect.IO
import cats.syntax.traverse._
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.{ErrorOr, _}
import common.validation.Parse.Parse
import cromwell.core._
import cromwell.languages.util.LanguageFactoryUtil
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import wdl.draft2.model.{WdlNamespace, WdlNamespaceWithWorkflow}
import wom.core.{WorkflowJson, WorkflowOptionsJson, WorkflowSource}
import wom.graph.GraphNodePort.OutputPort
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomExecutableMakers._
import wom.executable.WomBundle
import wom.transforms.WomExecutableMaker.ops._
import wom.values.WomValue

import scala.concurrent.Future

class WdlDraft2LanguageFactory() extends LanguageFactory {
  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                    workflowOptions: WorkflowOptions,
                                    importLocalFilesystem: Boolean,
                                    workflowIdForLogging: WorkflowId): Parse[ValidatedWomNamespace] = {

    def checkTypes(namespace: WdlNamespaceWithWorkflow, inputs: Map[OutputPort, WomValue]): Checked[Unit] = {
      val allDeclarations = namespace.workflow.declarations ++ namespace.workflow.calls.flatMap(_.declarations)
      val list: List[Checked[Unit]] = inputs.map({ case (k, v) =>
        allDeclarations.find(_.fullyQualifiedName == k) match {
          case Some(decl) if decl.womType.coerceRawValue(v).isFailure =>
            s"Invalid right-side type of '$k'.  Expecting ${decl.womType.toDisplayString}, got ${v.womType.toDisplayString}".invalidNelCheck[Unit]
          case _ => ().validNelCheck
        }
      }).toList

      list.sequence[Checked, Unit].void
    }

    val baseResolvers: List[String => WorkflowSource] = if (importLocalFilesystem) {
      List(WdlNamespace.fileResolver, WdlNamespace.httpResolver)
    } else {
      List(WdlNamespace.httpResolver)
    }

    import common.validation.Validation._

    val wdlNamespaceValidation: ErrorOr[WdlNamespaceWithWorkflow] = source match {
      case w: WorkflowSourceFilesWithDependenciesZip =>
        for {
          importsDir <- LanguageFactoryUtil.validateImportsDirectory(w.importsZip)
          betterFilesImportsDir = better.files.File(importsDir.pathAsString)
          directoryResolver = WdlNamespace.directoryResolver(betterFilesImportsDir): String => WorkflowSource
          resolvers = directoryResolver +: baseResolvers
          wf <- WdlNamespaceWithWorkflow.load(w.workflowSource, resolvers).toErrorOr
          _ = importsDir.delete(swallowIOExceptions = true)
        } yield wf
      case w: WorkflowSourceFilesWithoutImports =>
        WdlNamespaceWithWorkflow.load(w.workflowSource, baseResolvers).toErrorOr
    }

    def evaluateImports(wdlNamespace: WdlNamespace): Map[String, String] = {
      // Descend the namespace looking for imports and construct `MetadataEvent`s for them.
      def collectImportEvents: Map[String, String] = {
        (wdlNamespace.allNamespacesRecursively flatMap { ns =>
          ns.importUri.toList collect {
            // Do not publish import events for URIs which correspond to literal strings as these are the top-level
            // submitted workflow.
            case uri if uri != WdlNamespace.WorkflowResourceString => uri -> ns.sourceString
          }
        }).toMap
      }

      collectImportEvents
    }

    val checked: Checked[ValidatedWomNamespace] = for {
      wdlNamespace <- wdlNamespaceValidation.toEither
      _ <- validateWorkflowNameLengths(wdlNamespace)
      importedUris = evaluateImports(wdlNamespace)
      womExecutable <- wdlNamespace.toWomExecutable(Option(source.inputsJson))
      validatedWomNamespaceBeforeMetadata <- LanguageFactoryUtil.validateWomNamespace(womExecutable)
      _ <- checkTypes(wdlNamespace, validatedWomNamespaceBeforeMetadata.womValueInputs)
    } yield validatedWomNamespaceBeforeMetadata.copy(importedFileContent = importedUris)

    fromEither[IO](checked)
  }

  private def validateWorkflowNameLengths(namespace: WdlNamespaceWithWorkflow): Checked[Unit] = {
    import common.validation.Checked._
    def allWorkflowNames(n: WdlNamespace): Seq[String] = n.workflows.map(_.unqualifiedName) ++ n.namespaces.flatMap(allWorkflowNames)
    val tooLong = allWorkflowNames(namespace).filter(_.length >= 100)
    if (tooLong.nonEmpty) {
      ("Workflow names must be shorter than 100 characters: " + tooLong.mkString(" ")).invalidNelCheck
    } else {
      ().validNelCheck
    }
  }

  override def getWomBundle(workflowSource: WorkflowSource, workflowOptionsJson: WorkflowOptionsJson, importResolvers: List[String => Future[Checked[WomBundle]]]): Checked[WomBundle] =
    "getWomBundle method not implemented in WDL draft 2".invalidNelCheck

  override def createExecutable(womBundle: WomBundle, inputs: WorkflowJson): Checked[ValidatedWomNamespace] =
    "createExecutable method not implemented in WDL draft 2".invalidNelCheck
}
