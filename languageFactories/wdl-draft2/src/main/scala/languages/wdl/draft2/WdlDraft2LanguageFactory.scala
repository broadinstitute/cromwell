package languages.wdl.draft2

import cats.instances.either._
import cats.instances.list._
import cats.syntax.functor._
import cats.data.EitherT.fromEither
import cats.effect.IO
import cats.syntax.traverse._
import common.Checked
import common.validation.Validation._
import common.validation.Checked._
import common.validation.ErrorOr.{ErrorOr, _}
import common.validation.Parse.Parse
import cromwell.core._
import cromwell.languages.util.ImportResolver.ImportResolver
import cromwell.languages.util.{ImportResolver, LanguageFactoryUtil}
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import wdl.draft2.model.{Draft2ImportResolver, WdlNamespace, WdlNamespaceWithWorkflow}
import wdl.shared.transforms.wdlom2wom.WdlSharedInputParsing
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomBundleMakers._
import wom.core.{WorkflowJson, WorkflowOptionsJson, WorkflowSource}
import wom.graph.GraphNodePort.OutputPort
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomExecutableMakers._
import wom.executable.WomBundle
import wom.expression.IoFunctionSet
import wom.transforms.WomExecutableMaker.ops._
import wom.transforms.WomBundleMaker.ops._
import wom.values.WomValue
import languages.wdl.draft2.WdlDraft2LanguageFactory._

class WdlDraft2LanguageFactory(override val config: Map[String, Any]) extends LanguageFactory {

  override val languageName: String = "WDL"
  override val languageVersionName: String = "draft-2"

  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                    workflowOptions: WorkflowOptions,
                                    importLocalFilesystem: Boolean,
                                    workflowIdForLogging: WorkflowId,
                                    ioFunctions: IoFunctionSet): Parse[ValidatedWomNamespace] = {

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

    lazy val baseResolvers: List[String => WorkflowSource] = if (importLocalFilesystem) {
      List(WdlNamespace.fileResolver, httpResolver)
    } else {
      List(httpResolver)
    }

    import common.validation.Validation._

    lazy val wdlNamespaceValidation: ErrorOr[WdlNamespaceWithWorkflow] = source match {
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
      _ <- standardConfig.enabledCheck
      wdlNamespace <- wdlNamespaceValidation.toEither
      _ <- validateWorkflowNameLengths(wdlNamespace)
      importedUris = evaluateImports(wdlNamespace)
      womExecutable <- wdlNamespace.toWomExecutable(Option(source.inputsJson), ioFunctions, standardConfig.strictValidation)
      validatedWomNamespaceBeforeMetadata <- LanguageFactoryUtil.validateWomNamespace(womExecutable, ioFunctions)
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

  override def getWomBundle(workflowSource: WorkflowSource, workflowOptionsJson: WorkflowOptionsJson, importResolvers: List[ImportResolver], languageFactories: List[LanguageFactory]): Checked[WomBundle] = {
    for {
      _ <- standardConfig.enabledCheck
      namespace <- WdlNamespace.loadUsingSource(workflowSource, None, Some(importResolvers map resolverConverter)).toChecked
      womBundle <- namespace.toWomBundle
    } yield womBundle
  }

  override def createExecutable(womBundle: WomBundle, inputs: WorkflowJson, ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace] = for {
    _ <- standardConfig.enabledCheck
    executable <- WdlSharedInputParsing.buildWomExecutable(womBundle, Option(inputs), ioFunctions, standardConfig.strictValidation)
    validatedNamespace <- LanguageFactoryUtil.validateWomNamespace(executable, ioFunctions)
  } yield validatedNamespace

  // Commentary: we'll set this as the default in the reference.conf, so most people will get WDL draft 2 if nothing else looks parsable.
  override def looksParsable(content: String): Boolean = false
}

object WdlDraft2LanguageFactory {
  private def resolverConverter(importResolver: ImportResolver): Draft2ImportResolver = str => importResolver.run(str) match {
    case Right(imported) => imported
    case Left(errors) => throw new RuntimeException(s"Bad import $str: ${errors.toList.mkString(System.lineSeparator)}")
  }

  val httpResolver = resolverConverter(ImportResolver.httpResolver)
  def httpResolverWithHeaders(headers: Map[String, String]) = resolverConverter(ImportResolver.httpResolverWithHeaders(headers))
}
