package cromwell.backend.google.pipelines.v2alpha1

import cromwell.backend.google.pipelines.common.PipelinesApiAttributes.LocalizationConfiguration
import cromwell.backend.google.pipelines.common.PipelinesApiFileInput
import cromwell.backend.google.pipelines.common.io.{DiskType, PipelinesApiWorkingDisk}
import cromwell.core.path.DefaultPathBuilder
import cromwell.filesystems.demo.dos.DemoDosPathBuilder
import org.scalatest.{FlatSpec, Matchers}

import scala.collection.JavaConverters._
import eu.timepit.refined.refineMV

class PipelinesConversionsSpec extends FlatSpec with Matchers {

  behavior of "PipelinesConversions"
  implicit val localizationConfiguration = LocalizationConfiguration(refineMV(1))

  it should "create a Demo DOS input parameter" in {
    val demoDosPathBuilder = new DemoDosPathBuilder
    val demoDosPath = demoDosPathBuilder.build("dos://dos.example.org/aaaabbbb-cccc-dddd-eeee-abcd0000dcba").get
    val containerRelativePath = DefaultPathBuilder.get("path/to/file.bai")
    val mount = PipelinesApiWorkingDisk(DiskType.LOCAL, 1)
    val input = PipelinesApiFileInput("example", demoDosPath, containerRelativePath, mount)
    val actions = PipelinesConversions.inputToParameter.toActions(input, Nil)
    actions.size should be(1)
    val action = actions.head

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

    val labels = action.get("labels").asInstanceOf[java.util.Map[_, _]]
    labels.keySet.asScala should contain theSameElementsAs List("tag", "inputName")
    labels.get("tag") should be("Localization")
    labels.get("inputName") should be("example")
  }

}
