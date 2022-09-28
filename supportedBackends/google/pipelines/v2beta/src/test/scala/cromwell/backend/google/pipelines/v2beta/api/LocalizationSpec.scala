package cromwell.backend.google.pipelines.v2beta.api

import common.assertion.CromwellTimeoutSpec
import cromwell.backend.google.pipelines.common.PipelinesApiJobPaths.DrsLocalizationManifestName
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.jdk.CollectionConverters._

class LocalizationSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "Localization"

  it should "create the right action to localize DRS files using a manifest" in {

    val manifestPathString = s"path/to/${DrsLocalizationManifestName}"
    val manifestPath = DefaultPathBuilder.get(manifestPathString)
    val tagKey = "tag"
    val tagLabel = "myLabel"

    val action = Localization.drsAction(manifestPath, Nil, Map(tagKey -> tagLabel), None)
    action.keySet.asScala should contain theSameElementsAs
      Set("commands", "environment", "imageUri", "labels", "mounts")

    action.get("commands") should be(a[java.util.List[_]])
    action.get("commands").asInstanceOf[java.util.List[_]] should contain theSameElementsAs List(
      "-m", manifestPathString
    )

    action.get("mounts") should be(a[java.util.List[_]])
    action.get("mounts").asInstanceOf[java.util.List[_]] should be (empty)

    action.get("imageUri") should be("somerepo/drs-downloader:tagged")

    val actionLabels = action.get("labels").asInstanceOf[java.util.Map[_, _]]
    actionLabels.keySet.asScala should contain theSameElementsAs List("tag")
    actionLabels.get(tagKey) should be(tagLabel)
  }

  it should "create the right action to localize DRS files using a manifest with requester pays" in {

    val manifestPathString = s"path/to/${DrsLocalizationManifestName}"
    val manifestPath = DefaultPathBuilder.get(manifestPathString)
    val tagKey = "tag"
    val tagLabel = "myLabel"
    val requesterPaysProjectId = "123"

    val action = Localization.drsAction(manifestPath, Nil, Map(tagKey -> tagLabel), Option(requesterPaysProjectId))
    action.keySet.asScala should contain theSameElementsAs
      Set("commands", "environment", "imageUri", "labels", "mounts")

    action.get("commands") should be(a[java.util.List[_]])
    action.get("commands").asInstanceOf[java.util.List[_]] should contain theSameElementsAs List(
      "-m", manifestPathString, "-r", requesterPaysProjectId
    )

    action.get("mounts") should be(a[java.util.List[_]])
    action.get("mounts").asInstanceOf[java.util.List[_]] should be (empty)

    action.get("imageUri") should be("somerepo/drs-downloader:tagged")

    val actionLabels = action.get("labels").asInstanceOf[java.util.Map[_, _]]
    actionLabels.keySet.asScala should contain theSameElementsAs List("tag")
    actionLabels.get(tagKey) should be(tagLabel)
  }
}
