package centaur.cwl
import com.typesafe.config.Config
import io.circe.optics.JsonPath
import io.circe.optics.JsonPath._
import io.circe.yaml.Printer.StringStyle
import io.circe.{Json, JsonObject, yaml}
import net.ceedubs.ficus.Ficus._

/**
  * Tools to pre-process the CWL workflows and inputs before feeding them to Cromwell so they can be executed on PAPI.
  */
class PAPIPreprocessor(config: Config) {
  // GCS directory where inputs for conformance tests are stored
  private val gcsPrefix = {
    val rawPrefix = config.as[String]("papi.default-input-gcs-prefix")
    if (rawPrefix.endsWith("/")) rawPrefix else rawPrefix + "/"
  }

  // Default docker pull image
  val DefaultDockerPull = "dockerPull" -> Json.fromString("ubuntu:latest")
  
  // Default docker image to be injected in a pre-existing requirements array
  private val DefaultDockerRequirement: Json = {
    Json.obj(
      "class" -> Json.fromString("DockerRequirement"),
      DefaultDockerPull
    )
  }

  // Requirements array with default docker requirement
  private val DefaultDockerRequirementList: Json = {
    Json.obj(
      "requirements" -> Json.arr(DefaultDockerRequirement)
    )
  }

  // Parse value, apply f to it, and print it back to String using the printer
  private def process(value: String, f: Json => Json, printer: Json => String) = {
    yaml.parser.parse(value) match {
      case Left(error) => throw new Exception(error.getMessage)
      case Right(json) => printer(f(json))
    }
  }

  // Process and print back as YAML
  private def processYaml(value: String)(f: Json => Json) =
    process(value, f, yaml.Printer.spaces2.copy(stringStyle = StringStyle.DoubleQuoted).pretty)

  // Process and print back as JSON
  private def processJson(value: String)(f: Json => Json) = process(value, f, io.circe.Printer.spaces2.pretty)

  // Prefix the string at "key" with the gcs prefix
  private def prefixLocationWithGcs(key: String): Json => Json = root.selectDynamic(key).string.modify(gcsPrefix + _)

  // Prefix "location" and "default"
  private val prefix = prefixLocationWithGcs("location").compose(prefixLocationWithGcs("path"))

  // Function to check if the given json has the provided key / value pair
  private def hasKeyValue(key: String, value: String): Json => Boolean = {
    root.selectDynamic(key).string.exist(_.equalsIgnoreCase(value))
  }

  // Return true if the given json object represents a File
  private def isFile(obj: JsonObject) = hasKeyValue("class", "File")(Json.fromJsonObject(obj))

  // Return true if the given json object represents a Directory
  private def isDirectory(obj: JsonObject) = hasKeyValue("class", "Directory")(Json.fromJsonObject(obj))

  // Prefix the location or path in the json object if it's a file or directory, otherwise recurse over its fields
  private def prefixObject(obj: JsonObject) = {
    // If the object is file or a directory, prefix it with the gcs prefix
    if (isFile(obj) || isDirectory(obj)) {
      prefix(Json.fromJsonObject(obj))
        // Even if it's a file it may have secondary files. So keep recursing on its fields
        .mapObject(_.mapValues(prefixFiles))
    }
    // Otherwise recursively process its fields
    else Json.fromJsonObject(obj.mapValues(prefixFiles))
  }

  // Fold over the json recursively and prefix all files
  private def prefixFiles(json: Json): Json = json.fold(
    jsonNull = json,
    jsonBoolean = _ => json,
    jsonNumber = _ => json,
    jsonString = _ => json,
    jsonObject = prefixObject,
    jsonArray = arr => Json.arr(arr.map(prefixFiles): _*)
  )

  /**
    * Pre-process input file by prefixing all files and directories with the gcs prefix
    */
  def preProcessInput(input: String): String = processJson(input)(prefixFiles)

  // Check if the given path (as an array or object) has a DockerRequirement element
  def hasDocker(jsonPath: JsonPath)(json: Json): Boolean = {
    val hasDockerInArray: Json => Boolean = jsonPath.arr.exist(_.exists(hasKeyValue("class", "DockerRequirement")))
    val hasDockerInObject: Json => Boolean = jsonPath.obj.exist(_.kleisli("DockerRequirement").nonEmpty)

    hasDockerInArray(json) || hasDockerInObject(json)
  }

  // Check if the given Json has a docker image in hints or requirements
  def hasDocker(json: Json): Boolean = hasDocker(root.hints)(json) || hasDocker(root.requirements)(json)

  // Add a default docker requirement to the workflow if it doesn't have one
  private def addDefaultDocker(workflow: Json) = if (!hasDocker(workflow)) {
    /*
      * deepMerge does not combine objects together but replaces keys which would overwrite existing requirements
      * so first check if there are requirements already and if so add our docker one.
      * Also turns out that the requirements section can be either an array or an object.
      * When it gets saladed the object is transformed to an array but because we deal with unsaladed cwl here
      * we have to handle both cases.
     */
    val requirementsAsArray = root.requirements.arr.modifyOption(_ :+ DefaultDockerRequirement)(workflow)
    val requirementsAsObject = root.requirements.obj.modifyOption(_.add("DockerRequirement", Json.obj(DefaultDockerPull)))(workflow)

    requirementsAsArray
      .orElse(requirementsAsObject)
      .getOrElse(workflow.deepMerge(DefaultDockerRequirementList))
  } else workflow

  /**
    * Pre-process the workflow by adding a default docker hint iff it doesn't have one
    */
  def preProcessWorkflow(workflow: String) = processYaml(workflow) { json =>
    def process = addDefaultDocker _
    
    // Some files contain a list of tools / workflows under the "$graph" field. In this case recursively add docker default to them
    root.$graph.arr.modifyOption(_.map(process))(json)
      // otherwise just process the file as a single workflow / tool
      .getOrElse(process(json))
  }
}
