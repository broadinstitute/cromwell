package centaur.cwl
import better.files.File
import com.typesafe.config.Config
import common.util.StringUtil._
import common.validation.IOChecked.IOChecked
import cwl.preprocessor.CwlPreProcessor
import io.circe.optics.JsonPath
import io.circe.optics.JsonPath._
import io.circe.yaml.Printer.StringStyle
import io.circe.{Json, yaml}
import net.ceedubs.ficus.Ficus._
import wom.util.YamlUtils

/**
  * Tools to pre-process the CWL workflows and inputs before feeding them to Cromwell so they can be executed on PAPI.
  */
class CloudPreprocessor(config: Config, prefixConfigPath: String) {
  val cwlPreProcessor = new CwlPreProcessor()
  
  // Cloud directory where inputs for conformance tests are stored
  private val cloudPrefix = config.as[String](prefixConfigPath).ensureSlashed

  // Default docker pull image
  val DefaultDockerPull = "dockerPull" -> Json.fromString("ubuntu:latest")
  
  // Default docker image to be injected in a pre-existing requirements array
  private val DefaultDockerHint: Json = {
    Json.obj(
      "class" -> Json.fromString("DockerRequirement"),
      DefaultDockerPull
    )
  }

  // hints array with default docker requirement
  private val DefaultDockerHintList: Json = {
    Json.obj(
      "hints" -> Json.arr(DefaultDockerHint)
    )
  }

  // Parse value, apply f to it, and print it back to String using the printer
  private def process(value: String, f: Json => Json, printer: Json => String) = {
    YamlUtils.parse(value) match {
      case Left(error) => throw new Exception(error.getMessage)
      case Right(json) => printer(f(json))
    }
  }

  // Process and print back as YAML
  private def processYaml(value: String)(f: Json => Json) =
    process(value, f, yaml.Printer.spaces2.copy(stringStyle = StringStyle.DoubleQuoted).pretty)

  // Prefix the string at "key" with the cloud prefix
  private def prefixLocation(value: String): String = {
    cloudPrefix + File(value.stripPrefix("file://")).name
  }

  // Function to check if the given json has the provided key / value pair
  private def hasKeyValue(key: String, value: String): Json => Boolean = {
    root.selectDynamic(key).string.exist(_.equalsIgnoreCase(value))
  }

  /**
    * Pre-process input file by prefixing all files and directories with the cloud prefix
    */
  def preProcessInput(input: String): IOChecked[String] = cwlPreProcessor.preProcessInputFiles(input, prefixLocation)

  // Check if the given path (as an array or object) has a DockerRequirement element
  def hasDocker(jsonPath: JsonPath)(json: Json): Boolean = {
    val hasDockerInArray: Json => Boolean = jsonPath.arr.exist(_.exists(hasKeyValue("class", "DockerRequirement")))
    val hasDockerInObject: Json => Boolean = jsonPath.obj.exist(_.kleisli("DockerRequirement").nonEmpty)

    hasDockerInArray(json) || hasDockerInObject(json)
  }

  // Check if the given Json has a docker image in hints or requirements
  def hasDocker(json: Json): Boolean = hasDocker(root.hints)(json) || hasDocker(root.requirements)(json)

  // Add a default docker hint to the workflow if it doesn't have one
  private val addDefaultDocker: Json => Json = workflow => if (!hasDocker(workflow)) {
    /*
      * deepMerge does not combine objects together but replaces keys which would overwrite existing hints
      * so first check if there are hints already and if so add our docker one.
      * Also turns out that the hints section can be either an array or an object.
      * When it gets saladed the object is transformed to an array but because we deal with unsaladed cwl here
      * we have to handle both cases.
     */
    val hintsAsArray = root.hints.arr.modifyOption(_ :+ DefaultDockerHint)(workflow)
    val hintsAsObject = root.hints.obj.modifyOption(_.add("DockerRequirement", Json.obj(DefaultDockerPull)))(workflow)

    hintsAsArray
      .orElse(hintsAsObject)
      .getOrElse(workflow.deepMerge(DefaultDockerHintList))
  } else workflow
  
  private val prefixDefaultFilesInCwl = CwlPreProcessor.mapFilesAndDirectories(prefixLocation) _

  /**
    * Pre-process the workflow by adding a default docker hint iff it doesn't have one
    */
  def preProcessWorkflow(workflow: String): String = processYaml(workflow)(addDefaultDocker.andThen(prefixDefaultFilesInCwl))
}
