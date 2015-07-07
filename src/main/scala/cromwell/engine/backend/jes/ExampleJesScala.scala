package cromwell.engine.backend.jes

import java.io.{FileInputStream, InputStreamReader}
import java.nio.file.Paths
import java.util
import java.util.UUID
import scala.collection.JavaConverters._
import com.google.api.client.extensions.java6.auth.oauth2.AuthorizationCodeInstalledApp
import com.google.api.client.googleapis.auth.oauth2.GoogleAuthorizationCodeFlow.Builder
import com.google.api.client.googleapis.auth.oauth2.GoogleClientSecrets
import com.google.api.client.googleapis.extensions.java6.auth.oauth2.GooglePromptReceiver
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.client.util.store.FileDataStoreFactory
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model._
import com.google.common.base.Joiner

// FIXME: This shouldn't be merged in!!!!!
object ExampleJesScala extends App {
  val RootGenomicsUrl = "https://staging-genomics.sandbox.googleapis.com"
  val DockerId = "gcr.io/broad-dsde-dev/ces-incipient:27"
  val ProjectId = "broad-dsde-dev"
  val SecretsFile = "/Users/jgentry/.google_client_secret.json"
  val GoogleUser = "jgentry"
  val GceLogsPath = "gs://jgentry-jes-dsde/"
  val Command = "echo starting; echo `date`; ls /work ; echo a ; echo `date`; java -jar /job/task/picard.jar ViewSam HEADER_ONLY=True I=/work/input.bam | tee /work/out.txt ; echo b ; echo `date` ; ls -l /work"
  //val Command = "exit -1"
  val DataStoreDir = Paths.get(System.getProperty("user.home"), ".jes-google-alpha")
  val JsonFactory = JacksonFactory.getDefaultInstance
  val PublicBam = "gs://broad-dsde-dev-public/NA12878/NA12878.subset.mapped.bam"
  val Scopes = Vector(
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/devstorage.full_control",
    "https://www.googleapis.com/auth/devstorage.read_write",
    "https://www.googleapis.com/auth/compute"
  )

  // Start what he had in main
  val id = create()
  println(s"Created pipeline with ID $id")
  val runName = run(id)
  println(s"Created pipeline run with name $runName")
  while(!check(runName)) Thread.sleep(5000)
  // End what he had in main

  def create(): String = {
    val cpr = new CreatePipelineRequest
    cpr.setProjectId(ProjectId)
    cpr.setDocker(makeDocker())
    cpr.setResources(makeResources())
    cpr.setName("No-name")

    val inFile = new Parameter()
    inFile.setName("input_bam")
    inFile.setValue("/work/input.bam")
    inFile.setType("REFERENCE")

    val outFile = new Parameter()
    outFile.setName("output")
    outFile.setValue("/work/out.txt")
    outFile.setType("REFERENCE")

    cpr.setParameters(Vector(inFile, outFile).asJava)

    pipelines().create(cpr).execute().getPipelineId
  }

  def run(id: String): String = {
    val rpr = new RunPipelineRequest
    rpr.setPipelineId(id)
    rpr.setProjectId(ProjectId)

    rpr.setResources(makeResources())

    val sa = new ServiceAccount()
    sa.setEmail("default")

    sa.setScopes(Scopes.asJava)
    rpr.setServiceAccount(sa)

    val execId = UUID.randomUUID().toString
    val logging = new Logging()
    logging.setGcsPath(s"$GceLogsPath$execId/system.log")
    rpr.setLogging(logging)

    val inputs = new util.HashMap[String, String]()
    inputs.put("input_bam", PublicBam)
    rpr.setInputs(inputs)

    val outputs = new util.HashMap[String, String]()
    outputs.put("output", s"$GceLogsPath$execId/output.txt")
    rpr.setOutputs(outputs)

    pipelines().run(rpr).execute().getName
  }

  def check(opName: String): Boolean = {
    println(s"Checking for operation " + opName)
    check(genomics().operations().get(opName).execute())
  }

  def check(op: Operation): Boolean = {
    if (op.getDone) {
      val jesStatus = Option(op.getError)
      printCheckStatus(jesStatus, op)
      val mapJoiner = Joiner.on(',').withKeyValueSeparator("=")
      println(s"Metadata is ${mapJoiner.join(op.getMetadata)}")
      true
    } else {
      false
    }
  }

  private def printCheckStatus(status: Option[Status], op: Operation): Unit = {
    status match {
      case Some(s) =>
        println(s"Error for ${op.getName}")
        println(s"Error Code: ${s.getCode}")
        println(s"Message: ${s.getMessage}")
        println(s"Details: ${s.getDetails}")
      case None =>
        println("Job was successful!")
        println(s"Created: ${op.getMetadata.get("created")}")
        println(s"Started: ${op.getMetadata.get("started")}")
        println(s"Finished: ${op.getMetadata.get("started")}")
    }
  }

  private def pipelines(): Genomics#Pipelines = {
    genomics().pipelines()
  }

  private def genomics(): Genomics = {
    val httpTransport = GoogleNetHttpTransport.newTrustedTransport()
    val secretStream = new InputStreamReader(new FileInputStream(SecretsFile))
    val clientSecrets = GoogleClientSecrets.load(JsonFactory, secretStream)

    val dataStoreFactory = new FileDataStoreFactory(DataStoreDir.toFile)
    val flow = new Builder(httpTransport,
      JsonFactory,
      clientSecrets,
      Scopes.asJava).setDataStoreFactory(dataStoreFactory).build()
    val credential = new AuthorizationCodeInstalledApp(flow, new GooglePromptReceiver()).authorize(GoogleUser)
    val builder = new Genomics.Builder(GoogleNetHttpTransport.newTrustedTransport(), JsonFactory, credential)
    builder.setApplicationName("no-application-name")
    builder.setRootUrl(RootGenomicsUrl)
    builder.build
  }

  private def makeResources(): Resources = {
    val r = new Resources()
    r.setCpu(1L)
    r.setRamGb(2.toDouble)
    r.setPreemptible(false)

    r.setZones(Vector("us-central1-a").asJava)

    val localDisk = new Disk()
    localDisk.setSizeGb(100L)
    localDisk.setType("LOCAL_SSD")
    localDisk.setName("local-disk")
    r.setDisks(Vector(localDisk).asJava)

    r
  }

  private def makeDocker(): DockerExecutor = {
    val docker = new DockerExecutor()
    docker.setImage(DockerId)
    docker.setCmd(s"/bin/bash -c '$Command'")
    docker
  }


}
