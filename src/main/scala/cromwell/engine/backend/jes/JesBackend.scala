package cromwell.engine.backend.jes

import java.io.{FileInputStream, InputStreamReader, File}
import java.net.{URI, URL}
import java.nio.file.Paths

import com.google.api.client.extensions.java6.auth.oauth2.AuthorizationCodeInstalledApp
import com.google.api.client.googleapis.auth.oauth2.GoogleAuthorizationCodeFlow
import com.google.api.client.googleapis.auth.oauth2.GoogleClientSecrets
import com.google.api.client.googleapis.extensions.java6.auth.oauth2.GooglePromptReceiver
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.client.util.store.FileDataStoreFactory
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.ServiceAccount
import com.typesafe.config.ConfigFactory
import com.typesafe.scalalogging.LazyLogging
import cromwell.binding.WdlExpression._
import cromwell.binding._
import cromwell.binding.values._
import cromwell.engine.EngineFunctions
import cromwell.engine.backend.Backend
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.db.DataAccess
import scala.collection.JavaConverters._
import scala.concurrent.{Future, ExecutionContext}
import scala.util.{Failure, Success, Try}
import JesBackend._

object JesBackend {
  // FIXME: Until Chris codes up the real ones. For now just braindeadedly store stdout and stderr so we can string it along
  class TemporaryEngineFunctions(gcsPath: String, stdout: String, stderr: String) extends EngineFunctions {
    override protected def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = ???
    override protected def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] = ???
    override protected def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = Success(WdlFile(stdout.stripPrefix(s"$gcsPath/")))
    override protected def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = Success(WdlFile(stderr.stripPrefix(s"$gcsPath/")))
  }

  // FIXME: WTF is this thing? It's for the connection boilerplate
  val Scopes = Vector( // FIXME: Should be in package object? I believe so
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/devstorage.full_control",
    "https://www.googleapis.com/auth/devstorage.read_write",
    "https://www.googleapis.com/auth/compute"
  )

  // NOTE: Used in connection boilerplate. All of that could probably be made cleaner w/ generators or something
  val JesServiceAccount = new ServiceAccount().setEmail("default").setScopes(JesBackend.Scopes.asJava)

  // Decoration around WorkflowDescriptor to generate bucket names and the like
  implicit class JesWorkflowDescriptor(val descriptor: WorkflowDescriptor) extends AnyVal {
    def bucket = s"$CromwellExecutionBucket/${descriptor.name}/${descriptor.id.toString}"
    def callDir(call: Call) = s"$bucket/call-${call.name}"
  }

  // NOTE: Used for connection boilerplate, could probably be made cleaner
  lazy val GenomicsService = buildGenomics

  private lazy val JesConf = ConfigFactory.load.getConfig("backend").getConfig("jes")
  lazy val GenomicsUrl = new URL(JesConf.getString("rootUrl"))
  lazy val GoogleSecrets = Paths.get(JesConf.getString("secretsFile"))
  lazy val GoogleUser = JesConf.getString("user")
  lazy val GoogleProject = JesConf.getString("project")
  lazy val GoogleApplicationName = JesConf.getString("applicationName")

  val CromwellExecutionBucket = s"gs://cromwell-dev/cromwell-executions"

  /*
    FIXME: At least for now the only files that can be used are stdout/stderr. However this leads to a problem
    where stdout.txt is input and output :) Redirect stdout/stderr to a different name, but it'll be localized back
    in GCS as stdout/stderr. Yes, it's hacky.
   */
  val LocalStdout = "job.stdout.txt"
  val LocalStderr = "job.stderr.txt"

  // NOTE: Connection boilerplate
  def buildGenomics: Genomics = {
    val jsonFactory = JacksonFactory.getDefaultInstance
    val httpTransport = GoogleNetHttpTransport.newTrustedTransport
    val secretStream = new InputStreamReader(new FileInputStream(GoogleSecrets.toFile))
    val clientSecrets = GoogleClientSecrets.load(jsonFactory, secretStream)
    // FIXME: The following shouldn't be hardcoded
    val dataStoreFactory = new FileDataStoreFactory(new File(System.getProperty("user.home"), ".jes-google-alpha"))
    val flow = new GoogleAuthorizationCodeFlow.Builder(httpTransport,
      jsonFactory,
      clientSecrets,
      Scopes.asJava).setDataStoreFactory(dataStoreFactory).build
    val credential = new AuthorizationCodeInstalledApp(flow, new GooglePromptReceiver).authorize(GoogleUser)
    new Genomics.Builder(httpTransport, jsonFactory, credential).setApplicationName(GoogleApplicationName).setRootUrl(GenomicsUrl.toString).build
  }
}

class JesBackend extends Backend with LazyLogging {
  // FIXME: As of right now this isn't needed, but it probably will be as things get less hacky
  override def adjustInputPaths(call: Call, inputs: CallInputs): CallInputs = inputs

  // FIXME: Pass through for now, but perhaps not as time goes on
  override def initializeForWorkflow(workflow: WorkflowDescriptor): HostInputs = {
    // FIXME: No need to copy GCS inputs for the workflow we should be able to direclty reference them
    workflow.actualInputs
  }

  override def executeCommand(commandLine: String,
                              workflowDescriptor: WorkflowDescriptor,
                              call: Call,
                              backendInputs: CallInputs,
                              scopedLookupFunction: ScopedLookupFunction):Try[Map[String, WdlValue]] = {
    val gcsPath = s"${workflowDescriptor.callDir(call)}"

    // FIXME: These are only used for that TemporaryEngineFunctions at the moment
    val stdout = s"$gcsPath/stdout.txt"
    val stderr = s"$gcsPath/stderr.txt"

    // FIXME: Not possible to currently get stdout/stderr so redirect everything and hope the WDL isn't doing that too
    val redirectedCommand = s"$commandLine > $LocalStdout  2> $LocalStderr"

    val status = Pipeline(redirectedCommand, workflowDescriptor, call, backendInputs, GoogleProject, GenomicsService).run.waitUntilComplete()

    // FIXME: We don't have the possiblyAutoConvertedValue from LocalBackend, I think we need that?

    // FIXME: This is *only* going to work w/ my hacked stdout/stderr above, see LocalBackend for more info
    val outputMappings = call.task.outputs.map { taskOutput =>
      val rawValue = taskOutput.expression.evaluate(
        scopedLookupFunction,
        new TemporaryEngineFunctions(gcsPath, stdout, stderr)
      )
      println(s"JesBackend setting ${taskOutput.name} to $rawValue")
      taskOutput.name -> rawValue
    }.toMap

    status match {
      case Run.Success(created, started, finished) =>
        // FIXME: DRY cochise, this is C/P from LocalBackend
        val taskOutputEvaluationFailures = outputMappings.filter {_._2.isFailure}
        if (taskOutputEvaluationFailures.isEmpty) {
          val unwrappedMap = outputMappings.collect { case (name, Success(wdlValue) ) => name -> wdlValue }.toMap
          Success(unwrappedMap)
        } else {
          val message = taskOutputEvaluationFailures.collect { case (name, Failure(e))  => s"$name: $e" }.mkString("\n")
          Failure(new Throwable(s"Workflow ${workflowDescriptor.id}: $message"))
        }
      case Run.Failed(created, started, finished, errorCode, errorMessage) => // FIXME: errorMessage looks like it's just "pipeline run failed"?
        Failure(new Throwable(s"Workflow ${workflowDescriptor.id}: errorCode $errorCode for command: $commandLine. Message: $errorMessage"))
    }

  }

  override def handleCallRestarts(restartableWorkflows: Seq[RestartableWorkflow],
                                  dataAccess: DataAccess)(implicit ec: ExecutionContext): Future[Any] = Future("FIXME")
}
