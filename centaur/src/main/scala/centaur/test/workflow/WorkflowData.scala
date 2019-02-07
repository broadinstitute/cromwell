package centaur.test.workflow

import java.util.concurrent.Executors

import better.files.File
import cats.data.Validated._
import cats.effect.{ContextShift, IO}
import com.google.cloud.storage.{BlobId, StorageOptions}
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import configs.Result.{Failure, Success}
import configs.syntax._
import cromwell.api.model.Label
import net.ceedubs.ficus.Ficus._
import org.http4s.client.{Client, JavaNetClientBuilder}
import spray.json._

import scala.concurrent.ExecutionContext

/**
  * Stores workflow data / references.
  *
  * For some test suites the inputs and options may dynamically rendered with secure variables. For the other test
  * suites that don't require these inputs / options, the respective Function1[Unit, String] will not be invoked to read
  * the non-existent content.
  */
case class WorkflowData(workflowContent: Option[String],
                        workflowUrl: Option[String],
                        workflowRoot: Option[String],
                        workflowType: Option[String],
                        workflowTypeVersion: Option[String],
                        inputs: Option[IO[String]],
                        options: Option[IO[String]],
                        labels: List[Label],
                        zippedImports: Option[File],
                        secondOptions: Option[IO[String]] = None,
                        thirdOptions: Option[IO[String]] = None)

object WorkflowData {
  val blockingEC = ExecutionContext.fromExecutorService(Executors.newFixedThreadPool(5))
  implicit val cs: ContextShift[IO] = IO.contextShift(blockingEC)
  val httpClient: Client[IO] = JavaNetClientBuilder[IO](blockingEC).create
  val gcsStorage = StorageOptions.getDefaultInstance.getService

  def fromConfig(filesConfig: Config, fullConfig: Config, basePath: File): ErrorOr[WorkflowData] = {
    val workflowUrl = filesConfig.as[Option[String]]("workflowUrl")
    val workflowSourcePath = filesConfig.as[Option[String]]("workflow")

    (workflowSourcePath, workflowUrl) match {
      case (Some(workflowPath), None) => Valid(WorkflowData(
        workflowPath = Option(workflowPath),
        workflowUrl = None,
        filesConfig = filesConfig,
        fullConfig = fullConfig,
        basePath = basePath))
      case (None, Some(_)) => Valid(WorkflowData(
        workflowPath = None,
        workflowUrl = workflowUrl,
        filesConfig = filesConfig,
        fullConfig = fullConfig,
        basePath = basePath))
      case (Some(_), Some(_)) => invalidNel(s"Both 'workflow' path or 'workflowUrl' can't be provided.")
      case (None, None) => invalidNel(s"No 'workflow' path or 'workflowUrl' provided.")
    }
  }

  def apply(workflowPath: Option[String], workflowUrl: Option[String], filesConfig: Config, fullConfig: Config, basePath: File): WorkflowData = {
    def slurp(file: String): IO[String] = file match {
      case http if http.startsWith("http://") || http.startsWith("https://") => 
        httpClient.expect[String](http)
      case gcs if gcs.startsWith("gs://") =>
        val noScheme = gcs.stripPrefix("gs://")
        val firstSlashPosition = noScheme.indexOf("/")
        val blob = BlobId.of(noScheme.substring(0, firstSlashPosition), noScheme.substring(firstSlashPosition + 1))
        IO { gcsStorage.readAllBytes(blob).map(_.toChar).mkString }
      case local =>
        IO { basePath./(local).contentAsString }
    }
    
    def getOptionalFileContent(name: String): Option[IO[String]] = {
      filesConfig.getAs[String](name).map(slurp)
    }

    def getImports = filesConfig.get[List[String]]("imports") match {
      case Success(paths) => zipImports(paths map basePath./)
      case Failure(_) => None
    }

    def getImportsDirName(workflowPath: Option[File], workflowUrl: Option[String]): String = {
      workflowPath match {
        case Some(file) => file.name.replaceAll("\\.[^.]*$", "")
        case None => // workflow url is defined
          val fileName = workflowUrl.get.split("/").last
          fileName.replaceAll("\\.[^.]*$", "")
      }
    }

    def zipImports(imports: List[File]): Option[File] = {
      imports match {
        case Nil => None
        case _ =>
          val importsDirName = getImportsDirName(workflowPath.map(basePath./), workflowUrl)
          val importsDir = File.newTemporaryDirectory(importsDirName + "_imports")
          imports foreach { p =>
            val srcFile = File(p.pathAsString)
            val destFile = importsDir / srcFile.name
            srcFile.copyTo(destFile, overwrite = true)
          }

          Option(importsDir.zip())
      }
    }

    def getLabels: List[Label] = {
      import cromwell.api.model.LabelsJsonFormatter._
      getOptionalFileContent("labels") map { _.unsafeRunSync().parseJson.convertTo[List[Label]] } getOrElse List.empty
    }

    val workflowSource = workflowPath.map(slurp).map(_.unsafeRunSync())
    val workflowType = fullConfig.get[Option[String]]("workflowType").value
    val workflowTypeVersion = fullConfig.get[Option[String]]("workflowTypeVersion").value
    val workflowRoot = fullConfig.get[Option[String]]("workflowRoot").value

    WorkflowData(
      workflowContent = workflowSource,
      workflowUrl = workflowUrl,
      workflowRoot = workflowRoot,
      workflowType = workflowType,
      workflowTypeVersion = workflowTypeVersion,
      inputs = getOptionalFileContent("inputs"),
      options = getOptionalFileContent("options"),
      secondOptions = getOptionalFileContent(name = "second-options").orElse(getOptionalFileContent("options")),
      thirdOptions = getOptionalFileContent(name = "third-options").orElse(getOptionalFileContent("options")),
      labels = getLabels,
      zippedImports = getImports
    )
  }
}
