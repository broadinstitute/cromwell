package cwl.preprocessor

import cats.effect.{ContextShift, IO}
import cats.instances.list._
import cats.syntax.parallel._
import common.validation.ErrorOr.ErrorOr
import common.validation.IOChecked._
import common.validation.Validation._
import cwl.preprocessor.CwlReference.EnhancedCwlId
import cwl.preprocessor.CwlPreProcessor._
import io.circe.Json
import io.circe.optics.JsonPath._
import cwl.preprocessor.CwlCanonicalizer._

/**
  * The real guts of the CWL pre-processor is taking a CWL reference and producing a single, self-contained JSON from it.
  */
private [preprocessor] class CwlCanonicalizer(saladFunction: SaladFunction)(implicit cs: ContextShift[IO]) {

  def getCanonicalCwl(reference: CwlReference,
                      namespacesJsonOption: Option[Json] = None,
                      schemasJsonOption: Option[Json] = None): IOChecked[Json] = {
    flattenCwlReferenceInner(
      reference,
      Map.empty,
      Map.empty,
      Set.empty,
      namespacesJsonOption,
      schemasJsonOption).map(_.processedJson)
  }

  /**
    * Flatten the cwl reference given already known processed references.
    */
  private def flattenCwlReferenceInner(cwlReference: CwlReference,
                                       unProcessedReferences: UnProcessedReferences,
                                       processedReferences: ProcessedReferences,
                                       breadCrumbs: Set[CwlReference],
                                       namespacesJsonOption: Option[Json],
                                       schemasJsonOption: Option[Json]): IOChecked[ProcessedJsonAndDependencies] = {
    /*
     * Salad and parse from a CWL reference into a Json object
     */
    def saladAndParse(ref: CwlReference): IOChecked[Json] =  for {
      saladed <- saladFunction(ref)
      saladedJson <- parseJson(saladed)
    } yield saladedJson

    for {
      // parse the file containing the reference
      parsed <- saladAndParse(cwlReference)
      // Get a Map[CwlReference, Json] from the parsed file. If the file is a JSON object and only contains one node, the map will only have 1 element
      newUnProcessedReferences = mapIdToContent(parsed).toMap
      // The reference json in the file
      referenceJson <- newUnProcessedReferences
        .collectFirst({ case (ref, json) if ref.pointerWithinFile == cwlReference.pointerWithinFile => json })
        .toIOChecked(s"Cannot find a tool or workflow with ID '${cwlReference.pointerWithinFile}' in file ${cwlReference.pathAsString}'s set: [${newUnProcessedReferences.keySet.mkString(", ")}]")
      // Process the reference json
      processed <- flattenJson(
        referenceJson,
        newUnProcessedReferences ++ unProcessedReferences,
        processedReferences,
        breadCrumbs + cwlReference,
        namespacesJsonOption,
        schemasJsonOption
      )
    } yield processed
  }

  /**
    * Given a Json representing a tool or workflow, flattens it and return the other processed references that were generated.
    *
    * NB: Flatten here means two things:
    *  - Find references within the CWL and convert them into 'local' links
    *  - Create a map of canonical links to JSON CWL content
    *
    * @param saladedJson           json to process
    * @param unProcessedReferences references that have been parsed and saladed (we have the json), but not flattened yet
    * @param processedReferences   references that are fully flattened
    * @param breadCrumbs           list of references that brought us here
    * @param namespacesJsonOption  Namespaces from the original json
    * @param schemasJsonOption     Schemas from the original json
    */
  private def flattenJson(saladedJson: Json,
                          unProcessedReferences: UnProcessedReferences,
                          processedReferences: ProcessedReferences,
                          breadCrumbs: Set[CwlReference],
                          namespacesJsonOption: Option[Json],
                          schemasJsonOption: Option[Json]): IOChecked[ProcessedJsonAndDependencies] = {
    /*
     * Given a reference from a step's run field, flattens it and return it
     * @param unProcessedReferences references that have been parsed and saladed (we have the json), but not flattened yet.
     * @param checkedProcessedReferences references that are fully processed
     * @param cwlReference reference being processed
     * @return a new ProcessedReferences Map including this cwlReference processed along with all the dependencies
     *         that might have been processed recursively.
     */
    def processCwlRunReference(checkedProcessedReferences: IOChecked[ProcessedReferences],
                               cwlReference: CwlReference): IOChecked[ProcessedReferences] = {
      def processReference(processedReferences: ProcessedReferences) = {

        val result: IOChecked[ProcessedJsonAndDependencies] = unProcessedReferences.get(cwlReference) match {
          case Some(unProcessedReferenceJson) =>
            // Found the json in the unprocessed map, no need to reparse the file, just flatten this json
            flattenJson(
              unProcessedReferenceJson,
              unProcessedReferences,
              processedReferences,
              breadCrumbs,
              namespacesJsonOption,
              schemasJsonOption
            )
          case None =>
            // This is the first time we're seeing this reference, we need to parse its file and flatten it
            flattenCwlReferenceInner(
              cwlReference,
              unProcessedReferences,
              processedReferences,
              breadCrumbs,
              namespacesJsonOption,
              schemasJsonOption
            )
        }

        result map {
          // Return everything we've got (the previously known "processedReferences" + our new processed reference + everything that was processed to get to it)
          case ProcessedJsonAndDependencies(processed, newReferences) => processedReferences ++ newReferences + (cwlReference -> processed)
        }
      }

      def processIfNeeded(processedReferences: ProcessedReferences): IOChecked[ProcessedReferences] = {
        // If the reference has already been processed, no need to do anything
        if (processedReferences.contains(cwlReference)) processedReferences.validIOChecked
        // If the reference is in the bread crumbs it means we circled back to it: fail the pre-processing
        else if (breadCrumbs.contains(cwlReference)) s"Found a circular dependency on $cwlReference".invalidIOChecked
        // Otherwise let's see if we already have the json for it or if we need to process the file
        else processReference(processedReferences)
      }

      for {
        processedReferences <- checkedProcessedReferences
        newReferences <- processIfNeeded(processedReferences)
      } yield newReferences
    }

    def addJsonKeyValue(originalJson: Json, key: String, valueOption: Option[Json]): Json = {
      valueOption match {
        case Some(value) => originalJson.mapObject(_.add(key, value))
        case None => originalJson
      }
    }

    val namespacesJson = addJsonKeyValue(saladedJson, JsonKeyNamespaces, namespacesJsonOption)
    val schemasJson = addJsonKeyValue(namespacesJson, JsonKeySchemas, schemasJsonOption)

    import cats.syntax.apply._

    // Take the processed runs and inject them in the json
    def inlineProcessedJsons(newKnownReferences: ProcessedReferences, inlinedRunWorkflows: Map[String, ProcessedJsonAndDependencies]) = {
      // Provide a function to swap the run reference with its json content
      val lookupFunction: Json => Json = {
        json: Json => {
          val fromRunReferenceMap = for {
            asString <- json.asString
            reference <- asString.asReference
            embeddedJson <- newKnownReferences.get(reference)
          } yield embeddedJson

          val fromInlinedWorkflow = for {
            asObject <- json.asObject
            id <- asObject.kleisli("id")
            idAsString <- id.asString
            embeddedJson <- inlinedRunWorkflows.get(idAsString)
          } yield embeddedJson.processedJson

          fromRunReferenceMap.orElse(fromInlinedWorkflow).getOrElse(json)
        }
      }

      val flattenedJson = root.steps.each.run.json.modify(lookupFunction)(schemasJson)

      ProcessedJsonAndDependencies(flattenedJson, newKnownReferences ++ inlinedRunWorkflows.values.flatMap(_.processedDependencies))
    }

    /*
     * Given a json, collects all "steps.run" values that are JSON Strings, and convert them to CwlReferences.
     * A saladed JSON is assumed.
     */
    def findRunReferences(json: Json): List[CwlReference] = {
      json.asArray match {
        case Some(cwls) => cwls.toList.flatMap(findRunReferences)
        case _ => root.steps.each.run.string.getAll(json).flatMap(_.asReference).distinct
      }
    }

    /*
     * Given a json, collects all "steps.run" values that are JSON Objects representing a workflow.
     * A saladed JSON is assumed.
     * @return a Map[String, Json], where the key is the cwl id of the workflow, and the value its content
     */
    def findRunInlinedWorkflows(json: Json): ErrorOr[Map[String, Json]] = {
      import cats.instances.list._
      import cats.syntax.traverse._

      json.asArray match {
        case Some(cwls) => cwls.toList
          .flatTraverse(findRunInlinedWorkflows(_).map(_.toList))
          .map(_.toMap)
        case _ =>
          // Look for all the "run" steps that are json objects
          root.steps.each.run.obj.getAll(json)
            .map(Json.fromJsonObject)
            // Only keep the workflows (CommandLineTools don't have steps so no need to process them)
            .filter(root.`class`.string.exist(_.equalsIgnoreCase("Workflow")))
            .traverse[ErrorOr, (String, Json)]( obj =>
              // Find the id of the workflow
              root.id.string.getOption(obj)
                .toErrorOr("Programmer error: Workflow did not contain an id. Make sure the cwl has been saladed")
                .map(_ -> obj)
            ).map(_.toMap)
      }
    }

    // Recursively process the run references (where run is a string pointing to another Workflow / Tool)
    // TODO: it would be nice to accumulate failures here somehow (while still folding and be able to re-use
    // successfully processed references, so I don't know if ErrorOr would work)
    val processedRunReferences: IOChecked[ProcessedReferences] = findRunReferences(schemasJson).foldLeft(processedReferences.validIOChecked)(processCwlRunReference)

    // Recursively process the inlined run workflows (where run is a Json object representing a workflow)
    val processedInlineReferences: IOChecked[Map[String, ProcessedJsonAndDependencies]] = (for {
      inlineWorkflowReferences <- findRunInlinedWorkflows(saladedJson).toIOChecked
      flattenedWorkflows <- inlineWorkflowReferences.toList.parTraverse[IOChecked, IOCheckedPar, (String, ProcessedJsonAndDependencies)]({
        case (id, value) => flattenJson(value, unProcessedReferences, processedReferences, breadCrumbs, namespacesJsonOption, schemasJsonOption).map(id -> _)
      })
    } yield flattenedWorkflows).map(_.toMap)

    // Replace the unprocessed runs with their processed value
    (processedRunReferences, processedInlineReferences).tupled.map(Function.tupled(inlineProcessedJsons))
  }
}

private object CwlCanonicalizer {
  /**
    * A Cwl json that has been processed (saladed and flattened), as well as its processed dependencies.
    */
  final case class ProcessedJsonAndDependencies(processedJson: Json, processedDependencies: ProcessedReferences)

  final type UnProcessedReferences = Map[CwlReference, Json]
  final type ProcessedReferences = Map[CwlReference, Json]
}
