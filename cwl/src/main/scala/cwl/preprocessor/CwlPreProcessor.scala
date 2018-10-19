package cwl.preprocessor

import better.files.{File => BFile}
import cats.data.{EitherT, NonEmptyList}
import cats.effect.IO
import cats.instances.list._
import cats.syntax.either._
import cats.syntax.traverse._
import common.validation.ErrorOr.ErrorOr
import common.validation.Parse
import common.validation.Parse._
import common.validation.Validation._
import cwl.command.ParentName
import cwl.ontology.Schema
import cwl.preprocessor.CwlPreProcessor._
import cwl.{CwlDecoder, FileAndId, FullyQualifiedName}
import io.circe.optics.JsonPath._
import io.circe.{Json, JsonNumber, JsonObject}
import mouse.all._
import org.slf4j.LoggerFactory

object CwlPreProcessor {
  private val Log = LoggerFactory.getLogger("CwlPreProcessor")
  private val LocalScheme = "file://"

  private [preprocessor] object CwlReference {
    def fromString(in: String) = {
      in.startsWith(LocalScheme).option {
        FullyQualifiedName.maybeApply(in)(ParentName.empty) match {
          case Some(FileAndId(file, _, _)) => CwlReference(BFile(file.stripFilePrefix), in)
          case _ => CwlReference(BFile(in.stripFilePrefix), in)
        }
      }
    }

    def apply(file: BFile, pointer: Option[String]) = {
      // prepends file:// to the absolute file path
      val prefixedFile = s"$LocalScheme${file.pathAsString}"
      val fullReference = pointer.map(p => s"$prefixedFile#$p").getOrElse(prefixedFile)
      new CwlReference(file, fullReference)
    }
  }

  /**
    * Saladed CWLs reference other local CWL "node" (workflow or tool) using a URI as follow:
    * file:///path/to/file/containing/node.cwl[#pointer_to_node]
    * #pointer_to_node to node is optional, and will specify which workflow or tool is being targeted in the file.
    *
    * e.g:
    *   {
    *     "class": "Workflow",
    *     "id": "file:///path/to/workflow/workflow.cwl",
    *     ...
    *     "steps": [
    *       {
    *         "run": "file:///path/to/workflow/multi_tools.cwl#my_tool",
    *         ...
    *       }
    *     ]  
    *   }
    *
    * This snippet contains 2 references, one that is the ID of this workflow, the other one is the run step pointing to "my_tool" in "/path/to/workflow/multi_tools.cwl"
    *
    * @param file: the file containing the referenced node. e.g: File(/path/to/file/containing/node.cwl)
    * @param fullReference: the full reference string as it is found in the saladed json. e.g: "file:///path/to/file/containing/node.cwl#pointer_to_node"
    */
  private [preprocessor] case class CwlReference(file: BFile, fullReference: String) {
    override def toString = fullReference
    val pointer: Option[String] = fullReference.split("#") match {
      case Array(_, p) => Option(p)
      case _ => None
    }
  }

  private [preprocessor] type BreadCrumb = List[CwlReference]
  private [preprocessor] type ProcessedReferences = Map[CwlReference, Json]
  private [preprocessor] type UnProcessedReferences = Map[CwlReference, Json]

  /**
    * A Cwl json that has been processed (saladed and flattened), as well as its processed dependencies.
    */
  private case class ProcessedJsonAndDependencies(processedJson: Json, processedDependencies: ProcessedReferences)

  val saladCwlFile: BFile => Parse[String] = { file =>
    Log.info(s"Pre-Processing ${file.pathAsString}")
    CwlDecoder.saladCwlFile(file)
  }

  implicit class PrintableJson(val json: Json) extends AnyVal {
    def printCompact = io.circe.Printer.noSpaces.pretty(json)
  }

  private [preprocessor] implicit class EnhancedCwlId(val id: String) extends AnyVal {
    def asReference: Option[CwlReference] = CwlReference.fromString(id)
    def stripFilePrefix = id.stripPrefix(LocalScheme)
  }

  private [preprocessor] implicit class EnhancedCwlFromJsonRun(val jsonObject: JsonObject) extends AnyVal {
    def id: Option[String] = jsonObject.kleisli("id").get.asString
  }

  def noLogging = new CwlPreProcessor(CwlDecoder.saladCwlFile)
}

/**
  * Class to create a standalone version of a CWL file.
  * @param saladFunction function that takes a file and produce a saladed version of the content
  */
class CwlPreProcessor(saladFunction: BFile => Parse[String] = saladCwlFile) {
  // Modify the string at "key" using the mappingFunction
  private def mapStringValue(key: String, mappingFunction: String => String): Json => Json = root.selectDynamic(key).string.modify(mappingFunction)

  // Map "location" and "default"
  private def prefix(mappingFunction: String => String): Json => Json = mapStringValue("location", mappingFunction).compose(mapStringValue("path", mappingFunction))

  // Function to check if the given json has the provided key / value pair
  private def hasKeyValue(key: String, value: String): Json => Boolean = {
    root.selectDynamic(key).string.exist(_.equalsIgnoreCase(value))
  }

  // Return true if the given json object represents a File
  private def isFile(obj: JsonObject) = hasKeyValue("class", "File")(Json.fromJsonObject(obj))

  // Return true if the given json object represents a Directory
  private def isDirectory(obj: JsonObject) = hasKeyValue("class", "Directory")(Json.fromJsonObject(obj))

  // Prefix the location or path in the json object if it's a file or directory, otherwise recurse over its fields
  private def prefixObject(mappingFunction: String => String)(obj: JsonObject): Json = {
    // If the object is file or a directory, prefix it with the gcs prefix
    if (isFile(obj) || isDirectory(obj)) {
      prefix(mappingFunction)(Json.fromJsonObject(obj))
        // Even if it's a file it may have secondary files. So keep recursing on its fields
        .mapObject(_.mapValues(mapFilesAndDirectories(mappingFunction)))
    }
    // Otherwise recursively process its fields
    else Json.fromJsonObject(obj.mapValues(mapFilesAndDirectories(mappingFunction)))
  }

  // Fold over the json recursively and prefix all files
  def mapFilesAndDirectories(mappingFunction: String => String)(json: Json): Json = json.fold(
    jsonNull = json,
    jsonBoolean = _ => json,
    jsonNumber = _ => json,
    jsonString = _ => json,
    jsonObject = prefixObject(mappingFunction),
    jsonArray = arr => Json.arr(arr.map(mapFilesAndDirectories(mappingFunction)): _*)
  )

  private def mapNumbers(json: Json): Json = {
    // Circumvent Circe's scientific format for numbers: convert to a JSON String without exponential notation.
    def nonScientificNumberFormatting(jsonNumber: JsonNumber): Json = {
      val conversions = Stream[JsonNumber => Option[Any]](
        _.toBigInt.map(_.longValue()),
        _.toBigDecimal.map(_.doubleValue()),
        Function.const(Option("null")))

      // The `get` is safe because `Option("null")` guarantees a match even if the other two Stream elements
      // do not satisfy the predicate.
      conversions.map(_.apply(jsonNumber)).find(_.isDefined).flatten.get.toString |> Json.fromString
    }

    json.fold(
      jsonNull = json,
      jsonBoolean = _ => json,
      jsonNumber = nonScientificNumberFormatting,
      jsonString = _ => json,
      jsonObject = _.mapValues(mapNumbers) |> Json.fromJsonObject,
      jsonArray = _.map(mapNumbers) |> Json.fromValues
    )
  }

  def preProcessInputFiles(inputContent: String, mappingFunction: String => String): Parse[String] = for {
    parsed <- parseYaml(inputContent)
    mapped = parsed |> mapFilesAndDirectories(mappingFunction) |> mapNumbers
  } yield mapped.printCompact

  private val JsonKeyNamespaces = s"$$namespaces"
  private val JsonKeySchemas = s"$$schemas"

  /**
    * Pre-process a CWL file and create a standalone, runnable (given proper inputs), inlined version of its content.
    *
    * The general idea is to work on CwlReferences, starting from the one coming to this function in the form of file and optional root.
    * The goal is to look at the steps in this workflow that point to other references, and recursively flatten them until we can replace the step with
    * its flat version.
    *
    * There are 3 pieces of information that are carried around during this process:
    *  1) ProcessedReferences: A Map[CwlReference, Json] of CwlReference for which we have the fully processed (saladed AND flattened) Json value.
    *
    *  2) UnProcessedReferences: A Map[CwlReference, Json] of CwlReference for which we have the saladed but NOT flattened Json value.
    *     This can happen because a file can contain multiple tools / workflows. When we salad / parse this file, we get (CwlReference, Json) pairs
    *     for all the workflow / tools in the file, but they are not flattened yet.
    *     We keep this to avoid having to re-salad / re-parse files unnecessarily.
    *
    *  3) BreadCrumb: A List[CwlReference] used to follow the trail of CwlReferences that we are processing as we recurse down.
    *     This is used to be able to detect circular dependencies (if the cwlReference being processed is in that list, then we have a circular dependency) .
    *
    */
  def preProcessCwlFile(file: BFile, cwlRoot: Option[String]): Parse[Json] = {
    val reference = CwlReference(file, cwlRoot)

    def absoluteSchemaPaths(json: Json): Json = {
      json mapArray {
        _ map absoluteSchemaPaths
      } mapString {
        Schema.getIriPath(reference.fullReference, _)
      }
    }

    def flatten(json: Json): Parse[ProcessedJsonAndDependencies] = {
      val namespacesJsonOption: Option[Json] = json.asObject.flatMap(_.kleisli(JsonKeyNamespaces))
      val schemasJsonOption: Option[Json] = json.asObject.flatMap(_.kleisli(JsonKeySchemas)).map(absoluteSchemaPaths)
      flattenCwlReference(
        CwlReference(file, cwlRoot),
        Map.empty,
        Map.empty,
        Set.empty,
        namespacesJsonOption,
        schemasJsonOption
      )
    }

    def flattenOrByPass(json: Json): EitherT[IO, NonEmptyList[String], Json] = {
      val fileContentReference = for {
        asObject <- json.asObject
        fileContentId <- asObject.kleisli("id")
        stringId <- fileContentId.asString
        fileContentReference <- CwlReference.fromString(stringId)
      } yield fileContentReference

      fileContentReference match {
        // This by passes the pre-processing if the file already has an id for which the file part doesn't match the path of the file
        // passed to this function, as this would indicate that it has already been saladed and pre-processed.
        case Some(contentRef) if !contentRef.file.equals(reference.file) => json.validParse
        case _ => flatten(json).map(_.processedJson)
      }
    }

    for {
      original <- parseYaml(file.contentAsString)
      flattened <- flattenOrByPass(original)
    } yield flattened
  }

  /**
    * Convenience method to get the processed workflow as a string.
    */
  def preProcessCwlFileToString(file: BFile, cwlRoot: Option[String]): Parse[String] = {
    preProcessCwlFile(file, cwlRoot).map(_.printCompact)
  }

  /**
    * Flatten the cwl reference given already known processed references.
    */
  private def flattenCwlReference(cwlReference: CwlReference,
                                  unProcessedReferences: UnProcessedReferences,
                                  processedReferences: ProcessedReferences,
                                  breadCrumbs: Set[CwlReference],
                                  namespacesJsonOption: Option[Json],
                                  schemasJsonOption: Option[Json]
                                 ): Parse[ProcessedJsonAndDependencies] = {
    for {
      // parse the file containing the reference
      parsed <- saladAndParse(cwlReference.file)
      // Get a Map[CwlReference, Json] from the parsed file. If the file is a JSON object and only contains one node, the map will only have 1 element 
      newUnProcessedReferences = mapIdToContent(parsed).toMap
      // The reference json in the file
      referenceJson <- newUnProcessedReferences
        .collectFirst({ case (ref, json) if ref.pointer == cwlReference.pointer => json })
        .toParse(s"Cannot find a tool or workflow with ID ${cwlReference.fullReference} in file ${cwlReference.file.pathAsString}")
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
    * Given a reference from a step's run field, flattens it and return it
    * @param unProcessedReferences references that have been parsed and saladed (we have the json), but not flattened yet.
    * @param checkedProcessedReferences references that are fully processed
    * @param cwlReference reference being processed
    * @return a new ProcessedReferences Map including this cwlReference processed along with all the dependencies
    *         that might have been processed recursively.
    */
  private def processCwlRunReference(unProcessedReferences: UnProcessedReferences,
                                     breadCrumbs: Set[CwlReference],
                                     namespacesJsonOption: Option[Json],
                                     schemasJsonOption: Option[Json])
                                    (checkedProcessedReferences: Parse[ProcessedReferences],
                                     cwlReference: CwlReference): Parse[ProcessedReferences] = {
    def processReference(processedReferences: ProcessedReferences) = {

      val result: Parse[ProcessedJsonAndDependencies] = unProcessedReferences.get(cwlReference) match {
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
          flattenCwlReference(
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

    def processIfNeeded(processedReferences: ProcessedReferences): Parse[ProcessedReferences] = {
      // If the reference has already been processed, no need to do anything
      if (processedReferences.contains(cwlReference)) processedReferences.validParse
      // If the reference is in the bread crumbs it means we circled back to it: fail the pre-processing
      else if (breadCrumbs.contains(cwlReference)) s"Found a circular dependency on $cwlReference".invalidParse
      // Otherwise let's see if we already have the json for it or if we need to process the file
      else processReference(processedReferences)
    }

    for {
      processedReferences <- checkedProcessedReferences
      newReferences <- processIfNeeded(processedReferences)
    } yield newReferences
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
                          schemasJsonOption: Option[Json]): Parse[ProcessedJsonAndDependencies] = {
    val foldFunction = processCwlRunReference(
      unProcessedReferences,
      breadCrumbs,
      namespacesJsonOption,
      schemasJsonOption
    ) _

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
          } yield embeddedJson.json

          val fromInlinedWorkflow = for {
            asObject <- json.asObject
            id <- asObject.kleisli("id")
            idAsString <- id.asString
            embeddedJson <- inlinedRunWorkflows.get(idAsString)
          } yield embeddedJson.processedJson

          fromRunReferenceMap.orElse(fromInlinedWorkflow).getOrElse(json)
        }
      }

      val flattenedJson = root.steps.each.run.json.modify(lookupFunction)(schemasJson.json)
      
      ProcessedJsonAndDependencies(flattenedJson, newKnownReferences ++ inlinedRunWorkflows.values.flatMap(_.processedDependencies))
    }
    
    // Recursively process the run references (where run is a string pointing to another Workflow / Tool)
    // TODO: it would be nice to accumulate failures here somehow (while still folding and be able to re-use
    // successfully processed references, so I don't know if ErrorOr would work)     
    val processedRunReferences: Parse[ProcessedReferences] = findRunReferences(schemasJson.json).foldLeft(processedReferences.validParse)(foldFunction)
    
    // Recursively process the inlined run workflows (where run is a Json object representing a workflow)
    val processedInlineReferences: Parse[Map[String, ProcessedJsonAndDependencies]] = (for {
      inlineWorkflowReferences <- Parse.errorOrParse[Map[String, Json]] { findRunInlinedWorkflows(saladedJson.json) }
      flattenedWorkflows <- inlineWorkflowReferences.toList.traverse({
        case (id, value) => flattenJson(value, unProcessedReferences, processedReferences, breadCrumbs, namespacesJsonOption, schemasJsonOption).map(id -> _)
      })
    } yield flattenedWorkflows).map(_.toMap)

    // Replace the unprocessed runs with their processed value
    (processedRunReferences, processedInlineReferences).tupled.map(Function.tupled(inlineProcessedJsons))
  }

  /**
    * Salad and parse a string to Json
    */
  private def saladAndParse(file: BFile): Parse[Json] = for {
    saladed <- saladFunction(file)
    saladedJson <- parseJson(saladed)
  } yield saladedJson

  private def parseJson(in: String): Parse[Json] = {
    Parse.checkedParse(io.circe.parser.parse(in).leftMap(error => NonEmptyList.one(error.message)))
  }

  private def parseYaml(in: String): Parse[Json] = {
    Parse.checkedParse(io.circe.yaml.parser.parse(in).leftMap(error => NonEmptyList.one(error.message)))
  }

  /**
    * Given a json, collects all "steps.run" values that are JSON Strings, and convert them to CwlReferences.
    * A saladed JSON is assumed.
    */
  private def findRunReferences(json: Json): List[CwlReference] = {
    json.asArray match {
      case Some(cwls) => cwls.toList.flatMap(findRunReferences)
      case _ => root.steps.each.run.string.getAll(json).flatMap(_.asReference).distinct
    }
  }

  /**
    * Given a json, collects all "steps.run" values that are JSON Objects representing a workflow.
    * A saladed JSON is assumed.
    * @return a Map[String, Json], where the key is the cwl id of the workflow, and the value its content
    */
  private def findRunInlinedWorkflows(json: Json): ErrorOr[Map[String, Json]] = {
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
          .traverse( obj =>
          // Find the id of the workflow
          root.id.string.getOption(obj)
            .toErrorOr("Programmer error: Workflow did not contain an id. Make sure the cwl has been saladed")
            .map(_ -> obj)
        ).map(_.toMap)
    }
  }

  /**
    * Given a json, collect all tools or workflows and map them with their reference id.
    * A saladed JSON is assumed.
    */
  private def mapIdToContent(json: Json): List[(CwlReference, Json)] = {
    json.asArray match {
      case Some(cwls) => cwls.toList.flatMap(mapIdToContent)
      case None => root.id.string.getOption(json).flatMap(_.asReference).map(_ -> json).toList
    }
  }
}
