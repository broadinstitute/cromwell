package centaur.cwl

import better.files.File
import com.typesafe.config.Config
import common.validation.Parse.{Parse, _}
import cromwell.core.path.Obsolete.Paths
import cromwell.core.path.{DefaultPathBuilderFactory, PathBuilderFactory}
import cromwell.filesystems.gcs.GcsPathBuilderFactory
import cromwell.filesystems.ftp.FtpPathBuilderFactory
import cwl.preprocessor.CwlPreProcessor

sealed trait CentaurCwlRunnerRunMode {
  /**
    * Returns a user friendly of this run mode.
    */
  def description: String

  /**
    * Returns a factory that can create a path builder for this run mode.
    */
  def pathBuilderFactory: PathBuilderFactory

  /**
    * Runs any preprocessing as needed on a workflow.
    *
    * For example, may prefix relative paths so that absolute URLs are used.
    */
  def preProcessWorkflow(workflow: String): String = workflow

  /**
    * Runs any preprocessing as needed on inputs.
    *
    * For example, may prefix relative paths so that absolute URLs are used.
    */
  def preProcessInput(path: File): Parse[String] = path.contentAsString.validParse
}

object CentaurCwlRunnerRunMode {
  def fromConfig(conf: Config): CentaurCwlRunnerRunMode = {
    conf.getString("mode") match {
      case "local" => LocalRunMode
      case "papi" => PapiRunMode(conf)
      case "tesk" => TeskRunMode(conf)
      case unknown => throw new UnsupportedOperationException(s"mode not recognized: $unknown")
    }
  }
}

case object LocalRunMode extends CentaurCwlRunnerRunMode {
  val cwlPreProcessor = new CwlPreProcessor()
  override lazy val description: String = "local"
  override lazy val pathBuilderFactory: PathBuilderFactory = DefaultPathBuilderFactory
  /*
    * If the file is a relative local path, resolve it against the path of the input json.
   */
  private def inputFilesMapper(inputJsonPath: File)(file: String) = {
    if (!Paths.get(file).isAbsolute) inputJsonPath.sibling(file).toString
    else file
  }

  override def preProcessInput(input: File): Parse[String] = cwlPreProcessor.preProcessInputFiles(input.contentAsString, inputFilesMapper(input))
}

case class PapiRunMode(conf: Config) extends CentaurCwlRunnerRunMode {
  private lazy val googleConfig = conf.getConfig("google")
  private lazy val preprocessor = new CloudPreprocessor(conf, "papi.default-input-gcs-prefix")

  override lazy val description: String = s"papi ${googleConfig.getString("auth")}"

  override lazy val pathBuilderFactory: PathBuilderFactory = {
    GcsPathBuilderFactory(conf, googleConfig)
  }

  override def preProcessWorkflow(workflow: String): String = preprocessor.preProcessWorkflow(workflow)

  override def preProcessInput(input: File): Parse[String] = preprocessor.preProcessInput(input.contentAsString)
}

case class TeskRunMode(conf: Config) extends CentaurCwlRunnerRunMode {
  private lazy val ftpConfig = conf.getConfig("ftp")
  private lazy val preprocessor = new CloudPreprocessor(conf, "tesk.default-input-ftp-prefix")

  override lazy val description: String = s"tesk"

  override lazy val pathBuilderFactory: PathBuilderFactory = {
    new FtpPathBuilderFactory(conf, ftpConfig)
  }

  override def preProcessWorkflow(workflow: String): String = preprocessor.preProcessWorkflow(workflow)

  override def preProcessInput(input: File): Parse[String] = preprocessor.preProcessInput(input.contentAsString)
}
