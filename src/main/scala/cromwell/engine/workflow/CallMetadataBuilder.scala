package cromwell.engine.workflow

import cromwell.binding._
import cromwell.engine.backend.{CallMetadata, StdoutStderr}
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.db.slick._
import cromwell.engine.ExecutionIndex._
import cromwell.engine.SymbolStoreEntry
import org.joda.time.DateTime

import scala.language.postfixOps

/**
 * Builds call metadata suitable for return as part of a workflow metadata request.
 */
object CallMetadataBuilder {

  // The various data passed into the `build` method has a very different shape from what the workflow/metadata
  // endpoint needs to serve.  As an intermediate step to producing data in the correct format, this class aggregates
  // data by `ExecutionDatabaseKey`, which is a useful prerequisite for building results in the final data format.
  private case class AssembledCallMetadata(key: ExecutionDatabaseKey,
                                           execution: Execution,
                                           inputs: Traversable[SymbolStoreEntry] = Seq.empty,
                                           outputs: Option[Traversable[SymbolStoreEntry]] = None,
                                           streams: Option[StdoutStderr] = None,
                                           backend: Option[String] = None,
                                           jobId: Option[String] = None,
                                           backendStatus: Option[String] = None)

  // Types used in interim steps of the construction of the call metadata.
  // Map from an `ExecutionDatabaseKey` to the interim `AssembledCallMetadata` format.
  private type ExecutionMap = Map[ExecutionDatabaseKey, AssembledCallMetadata]
  // A function that transforms from one `ExecutionMap` to another.
  private type ExecutionMapTransformer = ExecutionMap => ExecutionMap

  object BackendValues {
    def extract(job: Any): BackendValues = {
      job match {
        case ji: LocalJob => BackendValues("Local")
        case ji: JesJob => BackendValues("JES", jobId = Option(ji.jesId.toString), status = ji.jesStatus)
        case ji: SgeJob => BackendValues("SGE", jobId = Option(ji.sgeJobNumber.toString))
      }
    }
  }
  // Case class to homogenize the datatypes for each supported backend.
  case class BackendValues(name: String, jobId: Option[String] = None, status: Option[String] = None)

  /**
   * Function to build the "base" transformer, this needs to run first in the sequence of transformers since it builds
   * the initial `ExecutionMap` entries.
   */
  private def buildBaseTransformer(executions: Traversable[Execution], executionKeys: Traversable[ExecutionDatabaseKey]): ExecutionMapTransformer =
  // The input ExecutionMap is ignored in this transformer!
    _ =>
      (for {
        execution <- executions
        key = ExecutionDatabaseKey(execution.callFqn, execution.index.toIndex)
        if !execution.callFqn.isScatter && !key.isCollector(executionKeys)
      } yield key -> AssembledCallMetadata(key, execution)).toMap


  /**
   * Function to build a transformer that adds inputs data to the entries in the input `ExecutionMap`.
   */
  private def buildInputsTransformer(callInputs: Traversable[SymbolStoreEntry]): ExecutionMapTransformer =
    executionMap => {
      // Remove collector entries for the purposes of this endpoint.
      val inputsNoCollectors = filterCollectorSymbols(callInputs)
      val inputsByKey = inputsNoCollectors.groupBy(i => ExecutionDatabaseKey(i.scope, i.key.index))

      for {
        (inputKey, inputs) <- inputsByKey
        // Inputs always have index None, but execution entries may have real indexes.  Look for executions by
        // FQN only, then give them copies of the inputs with matching indexes.
        (executionKey, assembledMetadata) <- executionMap
        if executionKey.fqn == inputKey.fqn
        indexedInputs = inputs map { case SymbolStoreEntry(key, wdlType, maybeWdlValue) =>
          new SymbolStoreEntry(key.copy(index = executionKey.index), wdlType, maybeWdlValue)
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
      val outputsByKey = outputsNoCollectors.groupBy(o => ExecutionDatabaseKey(o.scope, o.key.index))
      for {
        (key, outputs) <- outputsByKey
        baseMetadata = executionMap.get(key).get
      } yield key -> baseMetadata.copy(outputs = Option(outputs))
    }

  /**
   * Function to build a transformer that adds job data to the entries in the input `ExecutionMap`.
   */
  private def buildJobDataTransformer(executionKeys: Traversable[ExecutionDatabaseKey], jobMap: Map[ExecutionDatabaseKey, Any]): ExecutionMapTransformer =
    executionMap => {
      for {
        (key, job) <- jobMap.toTraversable
        if !key.fqn.isScatter && !key.isCollector(executionKeys)
        baseMetadata = executionMap.get(key).get
        backend = BackendValues.extract(job)
      } yield key -> baseMetadata.copy(backend = Option(backend.name), jobId = backend.jobId, backendStatus = backend.status)
    }.toMap

  /**
   * Function to build a transformer that adds standard streams data to the entries in the input `ExecutionMap`.
   */
  private def buildStreamsTransformer(standardStreamsMap: Map[FullyQualifiedName, Seq[StdoutStderr]]): ExecutionMapTransformer =
    executionMap => {
      val databaseKeysWithNoneIndexes = executionMap.keys groupBy { _.fqn } filter {
        case (fqn, edks) => edks.size == 1 && edks.head.index.isEmpty
      }

      def indexForFqn(fqn: FullyQualifiedName, rawIndex: Int): ExecutionIndex = {
        if (databaseKeysWithNoneIndexes.contains(fqn)) None else rawIndex.toIndex
      }

      for {
        (fqn, seqOfStreams) <- standardStreamsMap.toTraversable
        (streams, rawIndex) <- seqOfStreams.zipWithIndex
        key = ExecutionDatabaseKey(fqn, indexForFqn(fqn, rawIndex))
        baseMetadata = executionMap.get(key).get
      } yield key -> baseMetadata.copy(streams = Option(streams))
    }.toMap

  /** Remove symbols corresponding to collectors. */
  private def filterCollectorSymbols(symbols: Traversable[SymbolStoreEntry]): Traversable[SymbolStoreEntry] = {
    // Squash duplicate ExecutionDatabaseKeys.
    val databaseKeys = symbols map { s => ExecutionDatabaseKey(s.scope, s.key.index) } toSet
    // Look for FQNs with more than one ExecutionDatabaseKey per FQN.
    val fqnsWithCollectors = databaseKeys.groupBy(_.fqn).collect({ case (fqn, keys) if keys.size > 1 => fqn }).toSet
    // We know which FQNs with None indexes correspond to collectors, so filter matching symbols.
    symbols.filterNot { s => s.key.index.isEmpty && fqnsWithCollectors.contains(s.scope) }
  }

  /**
   *  Construct the map of `FullyQualifiedName`s to `Seq[CallMetadata]` that contains all the call metadata available
   *  for the specified parameters.
   */
  def build(executions: Traversable[Execution],
            standardStreamsMap: Map[FullyQualifiedName, Seq[StdoutStderr]],
            callInputs: Traversable[SymbolStoreEntry],
            callOutputs: Traversable[SymbolStoreEntry],
            jobMap: Map[ExecutionDatabaseKey, Any]): Map[FullyQualifiedName, Seq[CallMetadata]] = {

    val executionKeys = executions map { x => ExecutionDatabaseKey(x.callFqn, x.index.toIndex) }

    // Define a sequence of ExecutionMap transformer functions.
    val executionMapTransformers = Seq(
      buildBaseTransformer(executions, executionKeys),
      buildInputsTransformer(callInputs),
      buildOutputsTransformer(callOutputs),
      buildJobDataTransformer(executionKeys, jobMap),
      buildStreamsTransformer(standardStreamsMap)
    )

    // Fold a zero ExecutionMap across this Seq of functions.
    val executionMap = executionMapTransformers.foldLeft(Map.empty: ExecutionMap) {
      case (map, transformer) => map ++ transformer(map) }

    def symbolToMapEntry(symbol: Symbol) = {
      val clob = symbol.wdlValue.get
      symbol.name -> clob.getSubString(1, clob.length().toInt)
    }

    // Convert from the convenience AssembledCallMetadata format to the CallMetadata format
    // that the endpoint needs to serve up.
    def constructCallMetadata(metadata: AssembledCallMetadata): CallMetadata = {
      val inputsMap = metadata.inputs map { entry => entry.key.name -> entry.wdlValue.get } toMap
      val outputsMap = for {
        outputs <- metadata.outputs.toSeq
        entry <- outputs
      } yield entry.key.name -> entry.wdlValue.get

      CallMetadata(
        inputs = inputsMap,
        executionStatus = metadata.execution.status,
        backend = metadata.backend,
        backendStatus = metadata.backendStatus,
        outputs = Option(outputsMap.toMap),
        start = metadata.execution.startDt map { new DateTime(_) },
        end = metadata.execution.endDt map { new DateTime(_) },
        jobId = metadata.jobId,
        returnCode = metadata.execution.rc,
        stdout = metadata.streams map { _.stdout },
        stderr = metadata.streams map { _.stderr })
    }

    // The CallMetadatas need to be grouped by FQN and sorted within an FQN by index.
    for {
      (key, fqnGroupedMetadatas) <- executionMap.values.groupBy(_.key.fqn)
      fqnGroupedAndSortedMetadatas = fqnGroupedMetadatas.toSeq.sortBy(_.key.index)
    } yield key.fqn -> (fqnGroupedAndSortedMetadatas map constructCallMetadata)
  }
}
