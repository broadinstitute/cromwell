package languages.wdl.draft3

import akka.actor.ActorRef
import cats.data.EitherT.fromEither
import cats.effect.IO
import cats.syntax.traverse._
import common.Checked
import common.validation.ErrorOr.{ErrorOr, _}
import common.validation.Parse.Parse
import cromwell.core.CromwellGraphNode.CromwellEnhancedOutputPort
import cromwell.core._
import cromwell.engine.workflow.lifecycle.materialization.LanguageFactory
import cromwell.services.metadata.MetadataService.{PutMetadataAction, womValueToMetadataEvents}
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import languages.util.LanguageFactoryUtil
import wdl.draft3.Draft3VersionSpecifics
import wdl.{WdlNamespace, WdlNamespaceWithWorkflow}
import wom.core.WorkflowSource
import wom.executable.ValidatedWomNamespace
import wom.expression.IoFunctionSet
import wom.graph.GraphNodePort.OutputPort
import wom.values.WomValue

class WdlDraft3LanguageFactory() extends LanguageFactory {
  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                    workflowOptions: WorkflowOptions,
                                    ioFunctions: IoFunctionSet,
                                    importLocalFilesystem: Boolean,
                                    serviceRegistryActor: ActorRef,
                                    workflowIdForLogging: WorkflowId): Parse[ValidatedWomNamespace] = {
    import cats.instances.either._
    import cats.instances.list._
    import cats.syntax.either._
    import cats.syntax.functor._
    import common.validation.Checked._


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
          wf <- WdlNamespaceWithWorkflow.load(w.workflowSource, resolvers)(Draft3VersionSpecifics).toErrorOr
          _ = importsDir.delete(swallowIOExceptions = true)
        } yield wf
      case w: WorkflowSourceFilesWithoutImports =>
        WdlNamespaceWithWorkflow.load(w.workflowSource, baseResolvers)(Draft3VersionSpecifics).toErrorOr
    }

    /* Publish `MetadataEvent`s for all imports in this `WdlNamespace`. */
    def publishImportMetadata(wdlNamespace: WdlNamespace): Unit = {
      def metadataEventForImportedNamespace(ns: WdlNamespace): MetadataEvent = {
        import MetadataKey._
        import WorkflowMetadataKeys._
        // This should only be called on namespaces that are known to have a defined `importUri` so the .get is safe.
        val escapedUri = ns.importUri.get.escapeMeta
        MetadataEvent(MetadataKey(
          workflowIdForLogging, None, SubmissionSection, SubmissionSection_Imports, escapedUri), MetadataValue(ns.sourceString))
      }

      // Descend the namespace looking for imports and construct `MetadataEvent`s for them.
      def collectImportEvents: List[MetadataEvent] = {
        wdlNamespace.allNamespacesRecursively flatMap { ns =>
          ns.importUri.toList collect {
            // Do not publish import events for URIs which correspond to literal strings as these are the top-level
            // submitted workflow.
            case uri if uri != WdlNamespace.WorkflowResourceString => metadataEventForImportedNamespace(ns)
          }
        }
      }

      serviceRegistryActor ! PutMetadataAction(collectImportEvents)
    }

    val errorOr: ErrorOr[ValidatedWomNamespace] = (for {
      wdlNamespace <- wdlNamespaceValidation.toEither
      _ <- validateWorkflowNameLengths(wdlNamespace)
      _ = publishImportMetadata(wdlNamespace)
      womExecutable <- wdlNamespace.womExecutable(Option(source.inputsJson))
      validatedWomNamespace <- LanguageFactoryUtil.validateWomNamespace(womExecutable, ioFunctions)
      _ <- checkTypes(wdlNamespace, validatedWomNamespace.womValueInputs)
      _ = pushWfInputsToMetadataService(validatedWomNamespace.womValueInputs, serviceRegistryActor, workflowIdForLogging)
    } yield validatedWomNamespace).toValidated

    fromEither[IO](errorOr.toEither)
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

  private def pushWfInputsToMetadataService(workflowInputs: Map[OutputPort, WomValue],
                                            serviceRegistryActor: ActorRef,
                                            workflowIdForLogging: WorkflowId): Unit = {
    // Inputs
    val inputEvents = workflowInputs match {
      case empty if empty.isEmpty =>
        List(MetadataEvent.empty(MetadataKey(workflowIdForLogging, None,WorkflowMetadataKeys.Inputs)))
      case inputs =>
        inputs flatMap { case (outputPort, womValue) =>
          val inputName = outputPort.fullyQualifiedName
          womValueToMetadataEvents(MetadataKey(workflowIdForLogging, None, s"${WorkflowMetadataKeys.Inputs}:$inputName"), womValue)
        }
    }

    serviceRegistryActor ! PutMetadataAction(inputEvents)
  }

  // TODO WOM: resurect ?
  //  private def validateDeclarations(namespace: WdlNamespaceWithWorkflow,
  //                                   options: WorkflowOptions,
  //                                   coercedInputs: WorkflowCoercedInputs,
  //                                   pathBuilders: List[PathBuilder]): ErrorOr[WorkflowCoercedInputs] = {
  //    namespace.staticDeclarationsRecursive(coercedInputs, NoFunctions) match {
  //      case Success(d) => d.validNel
  //      case Failure(e) => s"Workflow has invalid declarations: ${e.getMessage}".invalidNel
  //    }
  //  }
}
