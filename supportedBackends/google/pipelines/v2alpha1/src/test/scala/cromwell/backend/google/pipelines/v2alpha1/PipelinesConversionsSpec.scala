package cromwell.backend.google.pipelines.v2alpha1

import java.nio.channels.ReadableByteChannel

import cats.effect.IO
import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, MarthaResponse}
import com.google.cloud.NoCredentials
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.LocalizationConfiguration
import cromwell.backend.google.pipelines.common.PipelinesApiFileInput
import cromwell.backend.google.pipelines.common.io.{DiskType, PipelinesApiWorkingDisk}
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder
import cromwell.core.path.DefaultPathBuilder
import cromwell.filesystems.drs.DrsPathBuilder
import eu.timepit.refined.refineMV
import org.apache.http.impl.client.HttpClientBuilder
import org.scalatest.{FlatSpec, Matchers}

import scala.collection.JavaConverters._

class PipelinesConversionsSpec extends FlatSpec with Matchers {

  behavior of "PipelinesConversions"
  implicit val localizationConfiguration = LocalizationConfiguration(refineMV(1))

  private val marthaConfig: Config = ConfigFactory.parseString(
    """martha {
      |   url = "http://matha-url"
      |   request.json-template = "{"key": "${holder}"}"
      |}
      |""".stripMargin
  )

  private lazy val fakeCredentials = NoCredentials.getInstance

  private lazy val httpClientBuilder = HttpClientBuilder.create()

  private def drsReadInterpreter(marthaResponse: MarthaResponse): IO[ReadableByteChannel] =
    throw new UnsupportedOperationException("Currently PipelinesConversionsSpec doesn't need to use drs read interpreter.")

  it should "create a DRS input parameter" in {

    val drsPathBuilder = DrsPathBuilder(
      new DrsCloudNioFileSystemProvider(marthaConfig, fakeCredentials, httpClientBuilder, drsReadInterpreter),
      None,
    )
    val drsPath = drsPathBuilder.build("dos://dos.example.org/aaaabbbb-cccc-dddd-eeee-abcd0000dcba").get
    val containerRelativePath = DefaultPathBuilder.get("path/to/file.bai")
    val mount = PipelinesApiWorkingDisk(DiskType.LOCAL, 1)
    val input = PipelinesApiFileInput("example", drsPath, containerRelativePath, mount)
    val actions = PipelinesConversions.inputToParameter.toActions(input, Nil)
    actions.size should be(2)

    val logging = actions.head

    logging.keySet.asScala should contain theSameElementsAs
      Set("commands", "flags", "imageUri", "labels", "mounts", "timeout")

    logging.get("commands") should be(a[java.util.List[_]])
    logging.get("commands").asInstanceOf[java.util.List[_]] should contain(
      """printf '%s %s\n' "$(date -u '+%Y/%m/%d %H:%M:%S')" """ +
        """Localizing\ input\ dos://dos.example.org/aaaabbbb-cccc-dddd-eeee-abcd0000dcba\ """ +
        """-\>\ /cromwell_root/path/to/file.bai"""
    )

    logging.get("flags") should be(a[java.util.List[_]])
    logging.get("flags").asInstanceOf[java.util.List[_]] should be (empty)

    logging.get("mounts") should be(a[java.util.List[_]])
    logging.get("mounts").asInstanceOf[java.util.List[_]] should be (empty)

    logging.get("imageUri") should be(ActionBuilder.cloudSdkImage)

    val loggingLabels = logging.get("labels").asInstanceOf[java.util.Map[_, _]]
    loggingLabels.keySet.asScala should contain theSameElementsAs List("logging", "inputName")
    loggingLabels.get("logging") should be("Localization")
    loggingLabels.get("inputName") should be("example")

    val action = actions.tail.head

    action.keySet.asScala should contain theSameElementsAs
      Set("commands", "environment", "imageUri", "labels", "mounts")

    action.get("commands") should be(a[java.util.List[_]])
    action.get("commands").asInstanceOf[java.util.List[_]] should contain theSameElementsAs List(
      "dos://dos.example.org/aaaabbbb-cccc-dddd-eeee-abcd0000dcba",
      "/cromwell_root/path/to/file.bai"
    )

    action.get("mounts") should be(a[java.util.List[_]])
    action.get("mounts").asInstanceOf[java.util.List[_]] should be (empty)

    action.get("imageUri") should be("somerepo/dos-downloader:tagged")

    val actionLabels = action.get("labels").asInstanceOf[java.util.Map[_, _]]
    actionLabels.keySet.asScala should contain theSameElementsAs List("tag", "inputName")
    actionLabels.get("tag") should be("Localization")
    actionLabels.get("inputName") should be("example")
  }

}
