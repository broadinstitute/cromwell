package cromwell.engine.workflow

import cromwell.database.obj.Execution
import cromwell.engine.ExecutionIndex._
import cromwell.engine._
import cromwell.engine.backend.jes.OldStyleJesBackend
import cromwell.engine.backend.local.OldStyleLocalBackend
import cromwell.engine.backend.sge.OldStyleSgeBackend
import cromwell.engine.backend.{CallLogs, OldStyleCallMetadata}
import cromwell.engine.db.EngineConverters.EnhancedExecution
import cromwell.engine.db.{ExecutionDatabaseKey, ExecutionInfosByExecution, ExecutionWithCacheData}
import org.joda.time.DateTime
import wdl4s._

import scala.language.postfixOps
import scala.util.Try

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class CallCacheData(allowResultReuse: Boolean, cacheHitWorkflow: Option[String], cacheHitCall: Option[String])

/**
 * Builds call metadata suitable for return as part of a workflow metadata request.
 */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object CallMetadataBuilder {

  // The various data passed into the `build` method has a very different shape from what the workflow/metadata
  // endpoint needs to serve.  As an intermediate step to producing data in the correct format, this class aggregates
  // data by `ExecutionDatabaseKey`, which is a useful prerequisite for building results in the final data format.
  private case class AssembledCallMetadata(key: ExecutionDatabaseKey,
                                           execution: Execution,
                                           inputs: Traversable[SymbolStoreEntry] = Seq.empty,
                                           outputs: Option[Traversable[SymbolStoreEntry]] = None,
                                           streams: Option[CallLogs] = None,
                                           backend: Option[String] = None,
                                           jobId: Option[String] = None,
                                           backendStatus: Option[String] = None,
                                           executionEvents: Seq[ExecutionEventEntry] = Seq.empty,
                                           runtimeAttributes: Map[String, String] = Map.empty,
                                           cache: Option[CallCacheData] = None,
                                           failures: Seq[FailureEventEntry] = Seq.empty)

  // Types used in interim steps of the construction of the call metadata.
  // Map from an `ExecutionDatabaseKey` to the interim `AssembledCallMetadata` format.
  private type ExecutionMap = Map[ExecutionDatabaseKey, AssembledCallMetadata]
  // A function that transforms from one `ExecutionMap` to another.
  private type ExecutionMapTransformer = ExecutionMap => ExecutionMap

  object BackendValues {
    // FIXME needs to traverse the list of pluggable backends rather than hardcoding a list.
    def extract(infos: ExecutionInfosByExecution): BackendValues = {
      def extractValue(key: String): Option[String] = infos.executionInfos collectFirst { case i if i.key == key => i.value } flatten

      infos.execution.backendType match {
        case "Local" => BackendValues("Local", jobId = extractValue(OldStyleLocalBackend.InfoKeys.Pid))
        case "JES" => BackendValues("JES", jobId = extractValue(OldStyleJesBackend.InfoKeys.JesRunId), status = extractValue(OldStyleJesBackend.InfoKeys.JesStatus))
        case "SGE" => BackendValues("SGE", jobId = extractValue(OldStyleSgeBackend.InfoKeys.JobNumber))
      }
    }
  }
  // Case class to homogenize the datatypes for each supported backend.
  case class BackendValues(name: String, jobId: Option[String] = None, status: Option[String] = None)

  private def findLastAttempt(executionMap: ExecutionMap, callFqn: String, index: ExecutionIndex) = {
    executionMap.keys.filter(k => k.fqn == callFqn && k.index == index).map(_.attempt).max
  }

  /**
   * Function to build the "base" transformer, this needs to run first in the sequence of transformers since it builds
   * the initial `ExecutionMap` entries.
   */
  private def buildBaseTransformer(executions: Traversable[Execution], executionKeys: Traversable[ExecutionDatabaseKey]): ExecutionMapTransformer =
  // The input ExecutionMap is ignored in this transformer!
    _ =>
      (for {
        execution <- executions
        key = execution.toKey
        if !key.isScatter && !key.isCollector(executionKeys)
      } yield key -> AssembledCallMetadata(key, execution)).toMap


  /**
   * Function to build a transformer that adds inputs data to the entries in the input `ExecutionMap`.
   */
  private def buildInputsTransformer(callInputs: Traversable[SymbolStoreEntry]): ExecutionMapTransformer =
    executionMap => {
      // Remove collector entries for the purposes of this endpoint.
      val inputsNoCollectors = filterCollectorSymbols(callInputs)

      for {
        (executionKey, assembledMetadata) <- executionMap
        // Inputs always have index None, but execution entries may have real indexes.  Look for executions by
        // FQN only, then give them copies of the inputs with matching indexes.
        indexedInputs = inputsNoCollectors.filter(_.scope == executionKey.fqn) map {
          case SymbolStoreEntry(key, wdlType, maybeWdlValue, maybeSymbolHash) =>
            new SymbolStoreEntry(key.copy(index = executionKey.index), wdlType, maybeWdlValue, maybeSymbolHash)
        }
      } yield executionKey -> assembledMetadata.copy(inputs = indexedInputs)
    }

  /**
   * Function to build a transformer that adds outputs data to the entries in the input `ExecutionMap`.
   */
  private def buildOutputsTransformer(callOutputs: Traversable[SymbolStoreEntry]): ExecutionMapTransformer =
    executionMap => {
      // Remove collector entries for the purposes of this endpoint.
      val outputsNoCollectors = filterCollectorSymbols(callOutputs)
      for {
        (executionKey, assembledMetadata) <- executionMap
        outputs = outputsNoCollectors filter { o =>
          val lastAttempt = findLastAttempt(executionMap, o.scope, o.key.index)
          ExecutionDatabaseKey(o.scope, o.key.index, lastAttempt) == executionKey
        }
      } yield executionKey -> assembledMetadata.copy(outputs = Option(outputs))
    }

  /**
    * Function to build a transformer that adds outputs data to the entries in the input `ExecutionMap`.
    */
  private def buildExecutionEventsTransformer(eventsMap: Map[ExecutionDatabaseKey, Seq[ExecutionEventEntry]]): ExecutionMapTransformer =
    executionMap => {
      for {
        (executionKey, assembledMetadata) <- executionMap
        events = eventsMap.getOrElse(executionKey, Seq.empty)
      } yield executionKey -> assembledMetadata.copy(executionEvents = events)
    }

  /**
    * Function to build a transformer that adds outputs data to the entries in the input `ExecutionMap`.
    */
  private def buildRuntimeAttributesTransformer(runtimeAttributesMap: Map[ExecutionDatabaseKey, Map[String, String]]): ExecutionMapTransformer =
    executionMap => {
      for {
        (executionKey, assembledMetadata) <- executionMap
        attrs = runtimeAttributesMap.getOrElse(executionKey, Map.empty)
      } yield executionKey -> assembledMetadata.copy(runtimeAttributes = attrs)
    }

  /**
    * Function to build a transformer that adds call cache data to the entries in the input `ExecutionMap`.
    */
  private def buildCallCacheTransformer(callCacheData: Traversable[ExecutionWithCacheData]): ExecutionMapTransformer =
    executionMap => {
      for {
        (executionKey, assembledMetadata) <- executionMap
        cacheData <- callCacheData
        if executionKey.fqn == cacheData.execution.callFqn
        metadata = CallCacheData(cacheData.execution.allowsResultReuse, cacheData.cacheHit.map(_.workflowId), cacheData.cacheHit.map(_.callName))
      } yield executionKey -> assembledMetadata.copy(cache = Option(metadata))
    }

  /**
    * Function to build a transformer that adds outputs data to the entries in the input `ExecutionMap`.
    */
  private def buildCallFailureTransformer(callFailureMap: Map[ExecutionDatabaseKey, Seq[FailureEventEntry]]): ExecutionMapTransformer =
    executionMap => {
      for {
        (executionKey, assembledMetadata) <- executionMap
        failureEvents = callFailureMap.getOrElse(executionKey, Seq.empty)
      } yield executionKey -> assembledMetadata.copy(failures = failureEvents)
    }

  /**
   * Function to build a transformer that adds job data to the entries in the input `ExecutionMap`.
   */
  private def buildExecutionInfoTransformer(executionsAndInfos: Traversable[ExecutionInfosByExecution],
                                            executionKeys: Traversable[ExecutionDatabaseKey]): ExecutionMapTransformer =
    executionMap => {
      for {
        (executionKey, assembledMetadata) <- executionMap
        if !executionKey.isScatter && !executionKey.isCollector(executionKeys)
        ei <- executionsAndInfos
        if ei.execution.toKey == executionKey
        e = ei.execution
        backendValues = BackendValues.extract(ei)
      } yield executionKey -> assembledMetadata.copy(
        backend = Option(e.backendType),
        jobId = backendValues.jobId,
        backendStatus = backendValues.status)
    }

  /**
   * Function to build a transformer that adds standard streams data to the entries in the input `ExecutionMap`.
   */
  private def buildStreamsTransformer(executionsAndInfos: Traversable[ExecutionInfosByExecution]):
  ExecutionMapTransformer = {
    executionMap => {
      for {
        (executionKey, assembledMetadata) <- executionMap
        callLogs = executionsAndInfos find {
          _.execution.toKey == executionKey
        } flatMap {
          _.callLogs
        }
      } yield executionKey -> assembledMetadata.copy(streams = callLogs)
    }
  }

  /** Remove symbols corresponding to collectors. */
  private def filterCollectorSymbols(symbols: Traversable[SymbolStoreEntry]): Traversable[SymbolStoreEntry] = {
    // Squash duplicate ExecutionDatabaseKeys.
    val databaseKeys = symbols map { s => ExecutionDatabaseKey(s.scope, s.key.index, 1) } toSet
    // Look for FQNs with at least one ExecutionStoreKey with defined index.
    val fqnsWithCollectors = databaseKeys.groupBy(_.fqn).collect({ case (fqn, keys) if keys.exists(_.index.isDefined) => fqn }).toSet
    // We know which FQNs with None indexes correspond to collectors, so filter matching symbols.
    symbols.filterNot { s => s.key.index.isEmpty && fqnsWithCollectors.contains(s.scope) }
  }

  /**
   *  Construct the map of `FullyQualifiedName`s to `Seq[CallMetadata]` that contains all the call metadata available
   *  for the specified parameters.
   */
  def build(infosByExecution: Traversable[ExecutionInfosByExecution],
            callInputs: Traversable[SymbolStoreEntry],
            callOutputs: Traversable[SymbolStoreEntry],
            executionEvents: Map[ExecutionDatabaseKey, Seq[ExecutionEventEntry]],
            runtimeAttributes: Map[ExecutionDatabaseKey, Map[String, String]],
            cacheData: Traversable[ExecutionWithCacheData],
            callFailures: Map[ExecutionDatabaseKey, Seq[FailureEventEntry]]): Map[FullyQualifiedName, Seq[OldStyleCallMetadata]] = {

    val executionKeys = infosByExecution map { _.execution.toKey }

    // Define a sequence of ExecutionMap transformer functions.
    val executionMapTransformers = Seq(
      buildBaseTransformer(infosByExecution map { _.execution }, executionKeys),
      buildInputsTransformer(callInputs),
      buildOutputsTransformer(callOutputs),
      buildExecutionEventsTransformer(executionEvents),
      buildExecutionInfoTransformer(infosByExecution, executionKeys),
      buildStreamsTransformer(infosByExecution),
      buildRuntimeAttributesTransformer(runtimeAttributes),
      buildCallFailureTransformer(callFailures),
      buildCallCacheTransformer(cacheData)
    )

    // Fold a zero ExecutionMap across this Seq of functions.
    val executionMap = executionMapTransformers.foldLeft(Map.empty: ExecutionMap) {
      case (map, transformer) => map ++ transformer(map) }

    // Convert from the convenience AssembledCallMetadata format to the CallMetadata format
    // that the endpoint needs to serve up.
    def constructCallMetadata(metadata: AssembledCallMetadata): OldStyleCallMetadata = {
      val inputsMap = metadata.inputs map { entry => entry.key.name -> entry.wdlValue.get } toMap
      val outputsMap = for {
        outputs <- metadata.outputs.toSeq
        entry <- outputs
      } yield entry.key.name -> entry.wdlValue.get

      val attempt = metadata.execution.attempt
      val preemptible = metadata.runtimeAttributes.get("preemptible") flatMap { x => Try(x.toInt).toOption } map { _ >= attempt }
      val failures = if (metadata.failures.isEmpty) None else Option(metadata.failures)

      OldStyleCallMetadata(
        inputs = inputsMap,
        executionStatus = metadata.execution.status,
        backend = metadata.backend,
        backendStatus = metadata.backendStatus,
        outputs = Option(outputsMap.toMap),
        start = metadata.execution.startDt map { new DateTime(_) },
        end = metadata.execution.endDt map { new DateTime(_) },
        jobId = metadata.jobId,
        returnCode = metadata.execution.rc,
        shardIndex = metadata.execution.index,
        stdout = metadata.streams map { _.stdout },
        stderr = metadata.streams map { _.stderr },
        backendLogs = metadata.streams flatMap { _.backendLogs },
        executionEvents = metadata.executionEvents,
        attempt = attempt,
        runtimeAttributes = metadata.runtimeAttributes,
        preemptible = preemptible,
        failures = failures,
        cache = metadata.cache
      )
    }

    // The CallMetadatas need to be grouped by FQN and sorted within an FQN by index and within an index by attempt.
    for {
      (key, fqnGroupedMetadatas) <- executionMap.values.groupBy(_.key.fqn)
      fqnGroupedAndSortedMetadatas = fqnGroupedMetadatas.toSeq.sortBy(md => (md.key.index, md.key.attempt))
    } yield key.fqn -> (fqnGroupedAndSortedMetadatas map constructCallMetadata)
  }
}
