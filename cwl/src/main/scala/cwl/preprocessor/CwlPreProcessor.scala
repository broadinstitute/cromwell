package cwl.preprocessor

import java.util.concurrent.Executors

import cats.data.NonEmptyList
import cats.effect.{ContextShift, IO}
import cats.syntax.either._
import common.validation.IOChecked._
import cwl.CwlDecoder
import cwl.ontology.Schema
import cwl.preprocessor.CwlPreProcessor._
import cwl.preprocessor.CwlReference.EnhancedCwlId
import io.circe.optics.JsonPath._
import io.circe.{Json, JsonNumber, JsonObject}
import mouse.all._
import org.slf4j.LoggerFactory
import wom.util.YamlUtils

import scala.concurrent.ExecutionContext

/**
  * Class to create a standalone version of a CWL file.
  *
  * NB: Want to use the pre-processor? Use preProcessCwl(ref: CwlReference)
  *
  * @param saladFunction function that takes a file and produce a saladed version of the content
  */
class CwlPreProcessor(saladFunction: SaladFunction = saladCwlFile) {

  private val ec: ExecutionContext = ExecutionContext.fromExecutor(Executors.newFixedThreadPool(5))
  private implicit val cs = IO.contextShift(ec)
  
  /**
    * This is THE main entry point into the CWL pre-processor. Takes a CWL reference and
    * returns a canonical JSON version with all references resolved.
    *
    * @param ref The reference to the CWL to pre-process
    * @return A canonical JSON representation of the CWL with all internal references expanded in-place
    */
  def preProcessCwl(ref: CwlReference): IOChecked[Json] = ref match {
    case file: CwlFileReference => preProcessCwlFile(file)
    case other => preProcessRemoteCwl(other)
  }

  /**
    * Convenience method to get the processed workflow as a string.
    */
  def preProcessCwlToString(cwlReference: CwlReference): IOChecked[String] = preProcessCwl(cwlReference).map(_.printCompact)

  def preProcessInputFiles(inputContent: String, mappingFunction: String => String): IOChecked[String] = for {
    parsed <- parseYaml(inputContent)
    mapped = parsed |> mapFilesAndDirectories(mappingFunction) |> mapNumbers
  } yield mapped.printCompact

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
    *  3) BreadCrumb: A Set[CwlReference] used to follow the trail of CwlReferences that we are processing as we recurse down.
    *     This is used to be able to detect circular dependencies (if the cwlReference being processed is in that set, then we have a circular dependency) .
    *
    */
  private def preProcessCwlFile(reference: CwlFileReference): IOChecked[Json] = {

    def absoluteSchemaPaths(json: Json): Json = {
      json mapArray {
        _ map absoluteSchemaPaths
      } mapString {
        Schema.getIriPath(reference.fullReference, _)
      }
    }

    // NB the JSON here is only used to decide whether or not to flatten. If we do decide to flatten we throw away the
    // json and request a canonical version from the CwlCanonicalizer.
    def flattenOrByPass(json: Json): IOChecked[Json] = {
      def flatten(json: Json): IOChecked[Json] = {
        val cwlReferenceFlattener = new CwlCanonicalizer(saladFunction)
        val namespacesJsonOption: Option[Json] = json.asObject.flatMap(_.kleisli(JsonKeyNamespaces))
        val schemasJsonOption: Option[Json] = json.asObject.flatMap(_.kleisli(JsonKeySchemas)).map(absoluteSchemaPaths)

        cwlReferenceFlattener.getCanonicalCwl(
          reference,
          namespacesJsonOption,
          schemasJsonOption
        )
      }

      def bypass(alreadyCanonicalJson: Json): IOChecked[Json] = alreadyCanonicalJson.validIOChecked

      val fileContentReference = for {
        asObject <- json.asObject
        fileContentId <- asObject.kleisli("id")
        stringId <- fileContentId.asString
        fileContentReference <- CwlReference.fromString(stringId)
      } yield fileContentReference

      fileContentReference match {
        // This by passes the pre-processing if the file already has an id for which the file part doesn't match the path of the file
        // passed to this function, as this would indicate that it has already been saladed and pre-processed.
        case Some(CwlFileReference(file, _)) if !file.equals(reference.file) => bypass(json)
        case _ => flatten(json)
      }
    }

    for {
      original <- parseYaml(reference.file.contentAsString)
      flattened <- flattenOrByPass(original)
    } yield flattened
  }

  // Like 'File', except that we don't read any contents before passing the path over to cwltool to canonicalize.
  private def preProcessRemoteCwl(reference: CwlReference)(implicit cs: ContextShift[IO]): IOChecked[Json] = {
    val cwlCanonicalizer = new CwlCanonicalizer(saladFunction)
    cwlCanonicalizer.getCanonicalCwl(reference)
  }
}

object CwlPreProcessor {
  private val Log = LoggerFactory.getLogger("CwlPreProcessor")

  private [preprocessor] type SaladFunction = CwlReference => IOChecked[String]
  private [preprocessor] val JsonKeyNamespaces = s"$$namespaces"
  private [preprocessor] val JsonKeySchemas = s"$$schemas"

  private def saladSpinner(doLogging: Boolean): SaladFunction = ref => {
    if (doLogging) {
      Log.info(s"Pre-Processing ${ref.pathAsString}")
    }

    CwlDecoder.saladCwlFile(ref)
  }

  private [preprocessor] val saladCwlFile: SaladFunction = saladSpinner(true)
  private val saladCwlFileWithoutLogging: SaladFunction = saladSpinner(false)

  implicit class PrintableJson(val json: Json) extends AnyVal {
    def printCompact = io.circe.Printer.noSpaces.pretty(json)
  }

  def noLogging = new CwlPreProcessor(saladCwlFileWithoutLogging)

  // Fold over a json recursively and prefix all files
  def mapFilesAndDirectories(mappingFunction: String => String)(json: Json): Json = {
    // Function to check if the given json has the provided key / value pair
    def hasKeyValue(key: String, value: String): Json => Boolean = {
      root.selectDynamic(key).string.exist(_.equalsIgnoreCase(value))
    }

    // Return true if the given json object represents a File
    def isFile(obj: JsonObject) = hasKeyValue("class", "File")(Json.fromJsonObject(obj))

    // Return true if the given json object represents a Directory
    def isDirectory(obj: JsonObject) = hasKeyValue("class", "Directory")(Json.fromJsonObject(obj))

    // Modify the string at "key" using the mappingFunction
    def mapStringValue(key: String, mappingFunction: String => String): Json => Json = root.selectDynamic(key).string.modify(mappingFunction)

    // Map "location" and "default"
    def prefix(mappingFunction: String => String): Json => Json = mapStringValue("location", mappingFunction).compose(mapStringValue("path", mappingFunction))

    // Prefix the location or path in the json object if it's a file or directory, otherwise recurse over its fields
    def prefixObject(mappingFunction: String => String)(obj: JsonObject): Json = {
      // If the object is file or a directory, prefix it with the gcs prefix
      if (isFile(obj) || isDirectory(obj)) {
        prefix(mappingFunction)(Json.fromJsonObject(obj))
          // Even if it's a file it may have secondary files. So keep recursing on its fields
          .mapObject(_.mapValues(mapFilesAndDirectories(mappingFunction)))
      }
      // Otherwise recursively process its fields
      else Json.fromJsonObject(obj.mapValues(mapFilesAndDirectories(mappingFunction)))
    }

    json.fold(
      jsonNull = json,
      jsonBoolean = _ => json,
      jsonNumber = _ => json,
      jsonString = _ => json,
      jsonObject = prefixObject(mappingFunction),
      jsonArray = arr => Json.arr(arr.map(mapFilesAndDirectories(mappingFunction)): _*)
    )
  }

  private [preprocessor] def mapNumbers(json: Json): Json = {
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

  private [preprocessor] def parseJson(in: String): IOChecked[Json] = {
    io.circe.parser.parse(in).leftMap(error => NonEmptyList.one(error.message)).toIOChecked
  }

  private [preprocessor] def parseYaml(in: String): IOChecked[Json] = {
    val yaml = YamlUtils.parse(in)
    yaml.leftMap(error => NonEmptyList.one(error.message)).toIOChecked
  }

  /**
    * Given a json, collect all tools or workflows and map them with their reference id.
    * A saladed JSON is assumed.
    */
  private [preprocessor] def mapIdToContent(json: Json): List[(CwlReference, Json)] = {
    json.asArray match {
      case Some(cwls) => cwls.toList.flatMap(mapIdToContent)
      case None => root.id.string.getOption(json).flatMap(_.asReference).map(_ -> json).toList
    }
  }
}
