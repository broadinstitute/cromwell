package cromwell.backend.google.batch.runnable

import common.assertion.CromwellTimeoutSpec
import cromwell.backend.google.batch.models.GcpBatchJobPaths
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.OptionValues._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.jdk.CollectionConverters._

class LocalizationSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "Localization"

  it should "create the right runnable to localize DRS files using a manifest" in {
    val manifestPathString = s"path/to/${GcpBatchJobPaths.DrsLocalizationManifestName}"
    val manifestPath = DefaultPathBuilder.get(manifestPathString)
    val tagKey = "tag"
    val tagLabel = "myLabel"

    val container = Localization.drsRunnable(manifestPath, Map(tagKey -> tagLabel), None).getContainer
    val fields = container.getAllFields.asScala
    fields.keySet.map(_.getName) should contain theSameElementsAs Set("commands", "image_uri")
    //      Set("commands", "environment", "imageUri", "labels", "mounts")

    fields
      .find(_._1.getName == "commands")
      .value
      .asInstanceOf[(_, java.util.List[_])]
      ._2 should contain theSameElementsAs List(
      "-m",
      manifestPathString
    )

    //    runnable.get("mounts") should be(a[java.util.List[_]])
    //    runnable.get("mounts").asInstanceOf[java.util.List[_]] should be (empty)
    //
    fields.find(_._1.getName == "image_uri").value.asInstanceOf[(_, String)]._2 should be(
      "somerepo/drs-downloader:tagged"
    )
    //
    //    val actionLabels = runnable.get("labels").asInstanceOf[java.util.Map[_, _]]
    //    actionLabels.keySet.asScala should contain theSameElementsAs List("tag")
    //    actionLabels.get(tagKey) should be(tagLabel)
  }

  it should "create the right runnable to localize DRS files using a manifest with requester pays" in {

    val manifestPathString = s"path/to/${GcpBatchJobPaths.DrsLocalizationManifestName}"
    val manifestPath = DefaultPathBuilder.get(manifestPathString)
    val tagKey = "tag"
    val tagLabel = "myLabel"
    val requesterPaysProjectId = "123"

    val container =
      Localization.drsRunnable(manifestPath, Map(tagKey -> tagLabel), Option(requesterPaysProjectId)).getContainer
    val fields = container.getAllFields.asScala
    fields.keySet.map(_.getName) should contain theSameElementsAs Set("commands", "image_uri")
    //      Set("commands", "environment", "imageUri", "labels", "mounts")

    fields
      .find(_._1.getName == "commands")
      .value
      .asInstanceOf[(_, java.util.List[_])]
      ._2 should contain theSameElementsAs List(
      "-m",
      manifestPathString,
      "-r",
      requesterPaysProjectId
    )

    //    runnable.get("mounts") should be(a[java.util.List[_]])
    //    runnable.get("mounts").asInstanceOf[java.util.List[_]] should be (empty)

    fields.find(_._1.getName == "image_uri").value.asInstanceOf[(_, String)]._2 should be(
      "somerepo/drs-downloader:tagged"
    )

    //    val actionLabels = runnable.get("labels").asInstanceOf[java.util.Map[_, _]]
    //    actionLabels.keySet.asScala should contain theSameElementsAs List("tag")
    //    actionLabels.get(tagKey) should be(tagLabel)
  }
}
