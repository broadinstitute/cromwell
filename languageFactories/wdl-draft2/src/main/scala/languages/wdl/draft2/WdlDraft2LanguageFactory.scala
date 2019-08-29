package languages.wdl.draft2

import java.util.concurrent.Callable

import cats.data.EitherT.fromEither
import cats.effect.IO
import cats.instances.either._
import cats.instances.list._
import cats.syntax.functor._
import cats.syntax.traverse._
import com.typesafe.config.Config
import com.typesafe.scalalogging.StrictLogging
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr._
import common.validation.IOChecked.IOChecked
import common.validation.Validation._
import cromwell.core._
import cromwell.languages.util.ImportResolver.{ImportResolutionRequest, ImportResolver}
import cromwell.languages.util.ParserCache.ParserCacheInputs
import cromwell.languages.util.{ImportResolver, LanguageFactoryUtil, ParserCache}
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import languages.wdl.draft2.WdlDraft2LanguageFactory._
import mouse.all._
import net.ceedubs.ficus.Ficus._
import wdl.draft2.Draft2ResolvedImportBundle
import wdl.draft2.model.{Draft2ImportResolver, WdlNamespace, WdlNamespaceWithWorkflow, WdlNamespaceWithoutWorkflow}
import wdl.shared.transforms.wdlom2wom.WdlSharedInputParsing
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomBundleMakers._
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomExecutableMakers._
import wom.ResolvedImportRecord
import wom.core.{WorkflowJson, WorkflowOptionsJson, WorkflowSource}
import wom.executable.WomBundle
import wom.expression.IoFunctionSet
import wom.graph.GraphNodePort.OutputPort
import wom.transforms.WomBundleMaker.ops._
import wom.transforms.WomExecutableMaker.ops._
import wom.values._

import scala.concurrent.duration._
import scala.language.postfixOps

class WdlDraft2LanguageFactory(override val config: Config) extends LanguageFactory with ParserCache[WdlNamespace] with StrictLogging {

  override val languageName: String = "WDL"
  override val languageVersionName: String = "draft-2"

  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                 workflowSource: WorkflowSource,
                                 workflowOptions: WorkflowOptions,
                                 importLocalFilesystem: Boolean,
                                 workflowIdForLogging: WorkflowId,
                                 ioFunctions: IoFunctionSet,
                                 importResolvers: List[ImportResolver]): IOChecked[ValidatedWomNamespace] = {

    def checkTypes(namespace: WdlNamespace, inputs: Map[OutputPort, WomValue]): Checked[Unit] = namespace match {

      case namespaceWithWorkflow: WdlNamespaceWithWorkflow =>
        val allDeclarations = namespaceWithWorkflow.workflow.declarations ++ namespaceWithWorkflow.workflow.calls.flatMap(_.declarations)
        val list: List[Checked[Unit]] = inputs.map({ case (k, v) =>
          allDeclarations.find(_.fullyQualifiedName == k) match {
            case Some(decl) if decl.womType.coerceRawValue(v).isFailure =>
              s"Invalid right-side type of '$k'.  Expecting ${decl.womType.stableName}, got ${v.womType.stableName}".invalidNelCheck[Unit]
            case _ => ().validNelCheck
          }
        }).toList

        list.sequence[Checked, Unit].void

      case _: WdlNamespaceWithoutWorkflow =>
        logger.error("Programmer Error: validateNamespace should never get called on WdlNamespaceWithoutWorkflow")
        "Cannot execute this WDL: no primary workflow provided".invalidNelCheck
    }

    def validationCallable = new Callable[ErrorOr[WdlNamespace]] {
      def call: ErrorOr[WdlNamespace] = WdlNamespaceWithWorkflow.load(workflowSource, importResolvers map resolverConverter).toErrorOr
    }

    lazy val wdlNamespaceValidation: ErrorOr[WdlNamespace] = retrieveOrCalculate(ParserCacheInputs(Option(workflowSource), None, None, importResolvers), validationCallable)

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
      _ <- enabledCheck
      wdlNamespace <- wdlNamespaceValidation.toEither
      _ <- validateWorkflowNameLengths(wdlNamespace)
      importedUris = evaluateImports(wdlNamespace)
      womExecutable <- wdlNamespace.toWomExecutable(Option(source.inputsJson), ioFunctions, strictValidation)
      validatedWomNamespaceBeforeMetadata <- LanguageFactoryUtil.validateWomNamespace(womExecutable, ioFunctions)
      _ <- checkTypes(wdlNamespace, validatedWomNamespaceBeforeMetadata.womValueInputs)
    } yield validatedWomNamespaceBeforeMetadata.copy(importedFileContent = importedUris)

    fromEither[IO](checked)
  }

  private def validateWorkflowNameLengths(namespace: WdlNamespace): Checked[Unit] = {
    import common.validation.Checked._
    def allWorkflowNames(n: WdlNamespace): Seq[String] = n.workflows.map(_.unqualifiedName) ++ n.namespaces.flatMap(allWorkflowNames)
    val tooLong = allWorkflowNames(namespace).filter(_.length >= 100)
    if (tooLong.nonEmpty) {
      ("Workflow names must be shorter than 100 characters: " + tooLong.mkString(" ")).invalidNelCheck
    } else {
      ().validNelCheck
    }
  }

  override def getWomBundle(workflowSource: WorkflowSource,
                            workflowSourceOrigin: Option[ResolvedImportRecord],
                            workflowOptionsJson: WorkflowOptionsJson,
                            importResolvers: List[ImportResolver],
                            languageFactories: List[LanguageFactory],
                            convertNestedScatterToSubworkflow : Boolean = true): Checked[WomBundle] = {
    lazy val validationCallable = new Callable[ErrorOr[WdlNamespace]] {
      def call: ErrorOr[WdlNamespace] = WdlNamespace.loadUsingSource(workflowSource, None, Some(importResolvers map resolverConverter)).toErrorOr
    }

    lazy val parserCacheInputs = ParserCacheInputs(Option(workflowSource), workflowSourceOrigin.map(_.importPath), None, importResolvers)

    for {
      _ <- enabledCheck
      namespace <- retrieveOrCalculate(parserCacheInputs, validationCallable).toEither
      womBundle <- namespace.toWomBundle
    } yield womBundle.copyResolvedImportRecord(womBundle, workflowSourceOrigin)
  }

  override def createExecutable(womBundle: WomBundle, inputs: WorkflowJson, ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace] = for {
    _ <- enabledCheck
    executable <- WdlSharedInputParsing.buildWomExecutable(womBundle, Option(inputs), ioFunctions, strictValidation)
    validatedNamespace <- LanguageFactoryUtil.validateWomNamespace(executable, ioFunctions)
  } yield validatedNamespace

  // Commentary: we'll set this as the default in the reference.conf, so most people will get WDL draft 2 if nothing else looks parsable.
  override def looksParsable(content: String): Boolean = false

  private[draft2] lazy val cacheConfig: Option[CacheConfig] = {
    // WDL version 2 namespace caching is now opt-in.
    for {
      _ <- enabled.option(())
      caching <- config.as[Option[Config]]("caching")
      cc <- CacheConfig.optionalConfig(caching, defaultConcurrency = 2, defaultSize = 1000L, defaultTtl = 20 minutes)
    } yield cc
  }
}

object WdlDraft2LanguageFactory {
  private def resolverConverter(importResolver: ImportResolver): Draft2ImportResolver = str => importResolver.resolver.run(ImportResolutionRequest(str, List.empty)) match {
    case Right(imported) => Draft2ResolvedImportBundle(imported.source, imported.resolvedImportRecord)
    case Left(errors) => throw new RuntimeException(s"Bad import $str: ${errors.toList.mkString(System.lineSeparator)}")
  }

  val httpResolver = resolverConverter(ImportResolver.HttpResolver())
  def httpResolverWithHeaders(headers: Map[String, String]) = resolverConverter(ImportResolver.HttpResolver(headers = headers))
}
