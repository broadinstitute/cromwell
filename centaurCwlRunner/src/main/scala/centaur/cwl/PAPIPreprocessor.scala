package centaur.cwl
import io.circe.optics.JsonPath
import io.circe.optics.JsonPath._
import io.circe.{Json, JsonObject, yaml}

/**
  * Tools to pre-process the CWL workflows and inputs before feeding them to Cromwell so they can be executed on PAPI.
  */
object PAPIPreprocessor {
  // GCS directory where inputs for conformance tests are stored
  private val gcsPrefix = "gs://centaur-cwl-conformance/cwl-inputs/"
  
  // Default docker image to inject if none is provided
  private val DefaultDocker: Json = {
    Json.obj(
       "class" -> Json.fromString("DockerRequirement"),
       "dockerPull" -> Json.fromString("ubuntu:latest")
    )
  }
  
  private val DefaultDockerRequirement: Json = {
    Json.obj(
      "requirements" -> Json.arr(DefaultDocker)
    )
  }
  
  // parse value, apply f to it, and print it back to String using the printer
  private def process(value: String, printer: Json => String, f: Json => Json) = {
    yaml.parser.parse(value) match {
      case Left(error) => throw new Exception(error.getMessage)
      case Right(json) => printer(f(json))
    }
  }

  // process and print back as YAML
  private def processYaml(value: String)(f: Json => Json) = process(value, yaml.Printer.spaces2.pretty, f)
  
  // process and print back as JSON
  private def processJson(value: String)(f: Json => Json) = process(value, io.circe.Printer.spaces2.pretty, f)

  // prefix the string at "key" with the gcs prefix
  private def prefixLocationWithGcs(key: String): Json => Json = root.selectDynamic(key).string.modify(gcsPrefix + _)
  
  // Prefix "location" and "path"
  private val prefix = prefixLocationWithGcs("location").compose(prefixLocationWithGcs("path"))
  
  // Function to check if the given json has the provided key / value pair
  private def hasKeyValue(key: String, value: String): Json => Boolean = {
    root.selectDynamic(key).string.exist(_.equalsIgnoreCase(value))
  }
  
  private def isFile(obj: JsonObject) = hasKeyValue("class", "File")(Json.fromJsonObject(obj))
  private def isDirectory(obj: JsonObject) = hasKeyValue("class", "Directory")(Json.fromJsonObject(obj))

  private def prefixObject(obj: JsonObject) = {
    // If the object is file a file or a directory, prefix it with the gcs prefix
    if (isFile(obj) || isDirectory(obj)) prefix(Json.fromJsonObject(obj))
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
    * Pre-process input file by prefixing all files with the gcs prefix
    */
  def preProcessInput(input: String): String = processJson(input)(prefixFiles)
  
  // Check if the given path (as an array) has an DockerRequirement element
  def hasDocker(jsonPath: JsonPath): Json => Boolean = jsonPath.arr.exist(_.exists(hasKeyValue("class", "DockerRequirement")))
  
  // Check if the given Json has a docker image in hints or requirements
  def hasDocker(json: Json): Boolean = hasDocker(root.hints)(json) || hasDocker(root.requirements)(json)
  
  // Add a default docker requirement to the workflow if it doesn't have one
  private def addDefaultDocker(workflow: Json) = if (!hasDocker(workflow)) {
    // deepMerge does not combine objects together but replaces keys, so first manually see if there are requirements
    // already and if so add our docker one
    root.requirements.arr.modifyOption(_ :+ DefaultDocker)(workflow)
      // otherwise add a new "requirements" field with our default docker
      .getOrElse(workflow.deepMerge(DefaultDockerRequirement))
  } else workflow

  /**
    * Pre-process the workflow by adding a default docker hint iff it doesn't have one
    */
  def preProcessWorkflow(workflow: String) = processYaml(workflow) { json =>
    // Some files contain a list of tools / workflows under the "$graph" field. In this case recursively add docker default to them
    root.$graph.arr.modifyOption(_.map(addDefaultDocker))(json)
      // otherwise just process the file as a single workflow / tool
      .getOrElse(addDefaultDocker(json))
  }
}
