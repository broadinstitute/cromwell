package cromwell.backend.google.pipelines.v2alpha1.api

import java.util

import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.v2alpha1.GenomicsFactory
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Labels.{Key, Value}
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}

import scala.collection.JavaConverters._

class ActionBuilderSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "ActionBuilder"

  val dockerRunActions = Table(
    ("description", "action", "command"),
    ("a cloud sdk action", ActionBuilder.cloudSdkAction, s"docker run ${GenomicsFactory.CloudSdkImage}"),
    ("a cloud sdk action with args",
      ActionBuilder.cloudSdkAction.setCommands(List("bash", "-c", "echo hello").asJava),
      s"docker run ${GenomicsFactory.CloudSdkImage} bash -c echo\\ hello"
    ),
    ("a cloud sdk action with quotes in the args",
      ActionBuilder.cloudSdkAction.setCommands(List("bash", "-c", "echo hello m'lord").asJava),
      s"docker run ${GenomicsFactory.CloudSdkImage} bash -c echo\\ hello\\ m\\'lord"
    ),
    ("a cloud sdk action with a newline in the args",
      ActionBuilder.cloudSdkAction.setCommands(List("bash", "-c", "echo hello\\\nworld").asJava),
      s"docker run ${GenomicsFactory.CloudSdkImage} bash -c echo\\ hello\\\\world"
    ),
    ("an action with multiple args",
      new Action()
        .setImageUri("ubuntu")
        .setEnvironment(Map("ENV" -> "dev").asJava)
        .setEntrypoint("")
        .setCommands(List("bash", "-c", "echo hello").asJava)
        .setFlags(List(ActionFlag.PublishExposedPorts, ActionFlag.AlwaysRun).map(_.toString).asJava)
        .setMounts(List(
          new Mount().setDisk("read-only-disk").setPath("/read/only/container").setReadOnly(true),
          new Mount().setDisk("read-write-disk").setPath("/read/write/container"),
        ).asJava)
        .setName("my_container_name")
        .setPidNamespace("host")
        .setPortMappings(Map("8008" -> Int.box(8000)).asJava),
      "docker run --name my_container_name" +
        " -v /mnt/read-only-disk:/read/only/container:ro -v /mnt/read-write-disk:/read/write/container" +
        " -e ENV:dev --pid=host -P -p 8008:8000 --entrypoint= ubuntu bash -c echo\\ hello"
    ),
  )

  forAll(dockerRunActions) { (description, action, command) =>
    it should s"convert $description" in {
      ActionBuilder.toDockerRun(action) should be(command)
    }
  }


  def memoryRetryExpectedCommand(lookupString: String): util.List[String] = {
    List(
      "/bin/sh",
      "-c",
      s"grep -E -q '$lookupString' /cromwell_root/stderr ; echo $$? > /cromwell_root/memory_retry_rc"
    ).asJava
  }

  val mounts = List(new Mount().setDisk("read-only-disk").setPath("/read/only/container"))
  val memoryRetryActionExpectedFlags = List(ActionFlag.AlwaysRun.toString).asJava
  val memoryRetryActionExpectedLabels = Map(Key.Tag -> Value.RetryWithMoreMemory).asJava

  it should "return cloud sdk action for one key in retry-with-double-memory" in {
    val lookupKeyList = List("OutOfMemory")
    val expectedCommand = memoryRetryExpectedCommand(lookupKeyList.mkString("|"))

    val action = ActionBuilder.checkForMemoryRetryAction(lookupKeyList, mounts)

    action.getCommands shouldBe expectedCommand
    action.getFlags shouldBe memoryRetryActionExpectedFlags
    action.getLabels shouldBe memoryRetryActionExpectedLabels
    action.getMounts shouldBe mounts.asJava
  }

  it should "return cloud sdk action for multiple keys in retry-with-double-memory" in {
    val lookupKeyList = List("OutOfMemory", "Killed", "Exit123")
    val expectedCommand = memoryRetryExpectedCommand(lookupKeyList.mkString("|"))

    val action = ActionBuilder.checkForMemoryRetryAction(lookupKeyList, mounts)

    action.getCommands shouldBe expectedCommand
    action.getFlags shouldBe memoryRetryActionExpectedFlags
    action.getLabels shouldBe memoryRetryActionExpectedLabels
    action.getMounts shouldBe mounts.asJava
  }
}
