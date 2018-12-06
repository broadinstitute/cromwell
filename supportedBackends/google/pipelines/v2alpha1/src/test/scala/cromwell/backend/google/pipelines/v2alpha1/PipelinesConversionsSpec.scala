package cromwell.backend.google.pipelines.v2alpha1

import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.google.pipelines.common.PipelinesApiAttributes.LocalizationConfiguration
import cromwell.backend.google.pipelines.common.PipelinesApiFileInput
import cromwell.backend.google.pipelines.common.io.{DiskType, PipelinesApiWorkingDisk}
import cromwell.core.path.DefaultPathBuilder
import cromwell.filesystems.drs.DrsPathBuilder
import eu.timepit.refined.refineMV
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

  it should "create a DRS input parameter" in {
    val drsPathBuilder = DrsPathBuilder(new DrsCloudNioFileSystemProvider(marthaConfig))
    val drsPath = drsPathBuilder.build("dos://dos.example.org/aaaabbbb-cccc-dddd-eeee-abcd0000dcba").get
    val containerRelativePath = DefaultPathBuilder.get("path/to/file.bai")
    val mount = PipelinesApiWorkingDisk(DiskType.LOCAL, 1)
    val input = PipelinesApiFileInput("example", drsPath, containerRelativePath, mount)
    val actions = PipelinesConversions.inputToParameter.toActions(input, Nil)
    actions.size should be(2)

    val logging = actions.head

    logging.keySet.asScala should contain theSameElementsAs
      Set("commands", "flags", "imageUri", "labels", "mounts")

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

    logging.get("imageUri") should be("google/cloud-sdk:slim")

    val loggingLabels = logging.get("labels").asInstanceOf[java.util.Map[_, _]]
    loggingLabels.keySet.asScala should contain theSameElementsAs List("logging", "inputName")
    loggingLabels.get("logging") should be("Localization")
    loggingLabels.get("inputName") should be("example")

    val action = actions.tail.head

    action.keySet.asScala should contain theSameElementsAs
      Set("commands", "entrypoint", "environment", "imageUri", "labels", "mounts")

    action.get("commands") should be(a[java.util.List[_]])
    action.get("commands").asInstanceOf[java.util.List[_]] should contain(
      "/path/to/some_executable " +
      "before args " +
      "dos://dos.example.org/aaaabbbb-cccc-dddd-eeee-abcd0000dcba " +
      "middle args " +
      "/cromwell_root/path/to/file.bai ends args"
    )

    action.get("entrypoint") should be("")

    action.get("mounts") should be(a[java.util.List[_]])
    action.get("mounts").asInstanceOf[java.util.List[_]] should be (empty)

    action.get("imageUri") should be("somerepo/dos-downloader:tagged")

    val actionLabels = action.get("labels").asInstanceOf[java.util.Map[_, _]]
    actionLabels.keySet.asScala should contain theSameElementsAs List("tag", "inputName")
    actionLabels.get("tag") should be("Localization")
    actionLabels.get("inputName") should be("example")
  }

}
