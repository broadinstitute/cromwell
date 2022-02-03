package cromwell.niochecksumtest

import akka.actor.{ActorRef, ActorSystem}
import akka.stream.ActorMaterializer
import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, GoogleDrsCredentials}
import com.typesafe.config.ConfigFactory
import cromwell.cloudsupport.gcp.auth.UserMode
import cromwell.core.WorkflowOptions
import cromwell.core.io.{AsyncIo, IoCommandBuilder, IoPromiseProxyActor}
import cromwell.engine.io.IoActor
import cromwell.engine.io.IoActor.IoConfig
import cromwell.filesystems.drs.{DrsPath, DrsReader}
import cromwell.services.ServiceRegistryActor
import net.ceedubs.ficus.Ficus.toFicusConfig
import spray.json.JsObject

import scala.concurrent.Await
import scala.concurrent.duration.DurationInt
import scala.util.{Failure, Success, Try}

/**
  * This class is set up to read a single (small!) DRS file into memory from our dev environment.
  * Ensure that the Google account you're using has permission to access the DRS file in the
  * environment you're targeting (use your test Google account in dev Terra!).
  * (1) Log into your Google account locally: `gcloud auth login [your username]`
  * (2) Log into Terra in your desired environment and make sure the repo you're using is linked
  * (3) Ensure that the Martha instance below is the right one.
  * (4) Set `drsPathToResolve` below to your DRS path
  * (5) Update `mySecretFile` below to point to the local adc.json file for your Google account
  * (6) "Run" NioChecksumTest using IntelliJ
  * (7) Stop NioChecksumTest using IntelliJ because it does not stop on its own :(
  */
object NioChecksumTest {

  /////////////////////
  // UPDATE THESE
  // V V V V V
  private val marthaUrl = "https://us-central1-broad-dsde-dev.cloudfunctions.net/martha_v3"
  private val drsPathToResolve = "drs://dg.4DFC:011a6a54-1bfe-4df9-ae24-990b12a812d3"
  private val mySecretFile = "/Users/[USER]/.config/gcloud/legacy_credentials/[GOOGLE_ACCOUNT_NAME]/adc.json"
  // ^ ^ ^ ^ ^
  /////////////////////

  private lazy val systemConfig = ConfigFactory.load()

  private lazy val drsConfig = ConfigFactory.parseString(
    s"""martha.url = "$marthaUrl"
      | access-token-acceptable-ttl = 1 hour
      |""".stripMargin
  )

  private lazy val authMode = UserMode("user auth", mySecretFile)
  private lazy val credentials = GoogleDrsCredentials(
    authMode.credentials(
      // From GoogleAccessTokenStrategy in CromwellDrsLocalizer
      List(
        "https://www.googleapis.com/auth/userinfo.email",
        "https://www.googleapis.com/auth/userinfo.profile"
      )
    ),
    drsConfig
  )

  def main(args: Array[String]): Unit = {
    implicit val context: ActorSystem = ActorSystem("NioChecksumTest")
    implicit val materializer: ActorMaterializer = ActorMaterializer()
    val serviceRegistryActor: ActorRef = context.actorOf(
      ServiceRegistryActor.props(systemConfig),
      name = "NioChecksumTestServiceRegistryActor"
    )

    val ioActor: ActorRef = context.actorOf(
      IoActor.props(
        ioConfig = systemConfig.as[IoConfig],
        serviceRegistryActor = serviceRegistryActor,
        applicationName = "NioChecksumTestIoActor"),
      name = "IoActor"
    )

    val ioPromiseProxyActor = context.actorOf(
      IoPromiseProxyActor.props(ioActor),
      name = "IoPromiseProxyActor"
    )

    try {
      println(resolveDrsFile(
        drsPathToResolve,
        new AsyncIo(ioPromiseProxyActor, IoCommandBuilder())
      ))
    } finally {
      context.stop(ioPromiseProxyActor)
      context.stop(ioActor)
      context.stop(serviceRegistryActor)
    }

    // Compiler gets mad if we don't use the result of the terminate call.
    println(s"Termination result is ${Await.result(context.terminate(), 60.seconds)}")
  }

  def resolveDrsFile(pathToResolve: String, asyncIo: AsyncIo): String = {

    val drsFilesystemProvider = new DrsCloudNioFileSystemProvider(
      drsConfig,
      credentials,
      DrsReader.readInterpreter(Option(authMode), WorkflowOptions(JsObject.empty), None)
    )

    val drsPath = DrsPath(drsFilesystemProvider.getCloudNioPath(pathToResolve), None)

    val fileContents = Try(
      Await.result(
        asyncIo.contentAsStringAsync(drsPath, Option(100000), failOnOverflow = false),
        60.seconds
      )
    )

    val message = fileContents match {
      case Success(contents) => s"Success! Downloaded a file of length ${contents.length}"
      case Failure(e) => s"Failed with error $e"
    }
    message
  }
}
