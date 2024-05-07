package cromwell.backend.google.batch.runnable

import com.google.cloud.batch.v1.{Runnable, Volume}
import common.assertion.CromwellTimeoutSpec
import cromwell.backend.google.batch.runnable.RunnableBuilder.EnhancedRunnableBuilder
import cromwell.backend.google.batch.runnable.RunnableLabels._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

import scala.jdk.CollectionConverters._

class RunnableBuilderSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "RunnableBuilder"

  private val dockerRunRunnables = Table(
    ("description", "runnable", "command"),
    ("a cloud sdk runnable", RunnableBuilder.cloudSdkRunnable, s"docker run ${RunnableUtils.CloudSdkImage}"),
    ("a cloud sdk runnable with args",
     RunnableBuilder.cloudSdkRunnable.withCommand("bash", "-c", "echo hello"),
     s"docker run ${RunnableUtils.CloudSdkImage} bash -c echo\\ hello"
    ),
    ("a cloud sdk runnable with quotes in the args",
     RunnableBuilder.cloudSdkRunnable.withCommand("bash", "-c", "echo hello m'lord"),
     s"docker run ${RunnableUtils.CloudSdkImage} bash -c echo\\ hello\\ m\\'lord"
    ),
    ("a cloud sdk runnable with a newline in the args",
     RunnableBuilder.cloudSdkRunnable.withCommand("bash", "-c", "echo hello\\\nworld"),
     s"docker run ${RunnableUtils.CloudSdkImage} bash -c echo\\ hello\\\\world"
    ),
    ("an runnable with multiple args",
     Runnable
       .newBuilder()
       .setContainer(Runnable.Container.newBuilder.setImageUri("ubuntu"))
       .withEntrypointCommand("")
       .withCommand("bash", "-c", "echo hello")
       .withAlwaysRun(true)
       .withVolumes(
         List(
           Volume
             .newBuilder()
             .setDeviceName("read-only-disk")
             .setMountPath("/mnt/read/only/container")
             .addMountOptions("ro"),
           Volume
             .newBuilder()
             .setDeviceName("read-write-disk")
             .setMountPath("/mnt/read/write/container")
             .addMountOptions("rw")
         ).map(_.build())
       ),
     "docker run" +
       " -v /mnt/read/only/container:/mnt/read/only/container -v /mnt/read/write/container:/mnt/read/write/container" +
       " ubuntu bash -c echo\\ hello"
    )
  )

  forAll(dockerRunRunnables) { (description, runnable, command) =>
    it should s"convert $description" in {
      RunnableBuilder.toDockerRun(runnable) should be(command)
    }
  }

  private val memoryRetryExpectedEntrypoint = "/bin/sh"

  def memoryRetryExpectedCommand(lookupString: String): List[String] =
    List(
      "-c",
      s"grep -E -q '$lookupString' /cromwell_root/stderr ; echo $$? > /cromwell_root/memory_retry_rc"
    )

  val volumes = List(
    Volume
      .newBuilder()
      .setDeviceName("read-only-disk")
      .setMountPath("/mnt/read/only/container") /*.addMountOptions("ro")*/
  ).map(_.build())

  private val memoryRetryRunnableExpectedLabels = Map(Key.Tag -> Value.RetryWithMoreMemory).asJava

  it should "return cloud sdk runnable for one key in retry-with-double-memory" in {
    val lookupKeyList = List("OutOfMemory")
    val expectedCommand = memoryRetryExpectedCommand(lookupKeyList.mkString("|"))

    val runnable = RunnableBuilder.checkForMemoryRetryRunnable(lookupKeyList, volumes)

    runnable.getContainer.getEntrypoint shouldBe memoryRetryExpectedEntrypoint
    runnable.getContainer.getCommandsList.asScala shouldBe expectedCommand
    runnable.getAlwaysRun shouldBe true
    runnable.getLabelsMap shouldBe memoryRetryRunnableExpectedLabels
    runnable.getContainer.getVolumesList.asScala.toList shouldBe volumes.map(v =>
      s"${v.getMountPath}:${v.getMountPath}"
    )
  }

  it should "return cloud sdk runnable for multiple keys in retry-with-double-memory" in {
    val lookupKeyList = List("OutOfMemory", "Killed", "Exit123")
    val expectedCommand = memoryRetryExpectedCommand(lookupKeyList.mkString("|"))

    val runnable = RunnableBuilder.checkForMemoryRetryRunnable(lookupKeyList, volumes)

    runnable.getContainer.getEntrypoint shouldBe memoryRetryExpectedEntrypoint
    runnable.getContainer.getCommandsList.asScala shouldBe expectedCommand
    runnable.getAlwaysRun shouldBe true
    runnable.getLabelsMap shouldBe memoryRetryRunnableExpectedLabels
    runnable.getContainer.getVolumesList.asScala.toList shouldBe volumes.map(v =>
      s"${v.getMountPath}:${v.getMountPath}"
    )
  }
}
