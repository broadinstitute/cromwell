package cromwell.backend.google.pipelines.v2beta

import java.nio.file.Paths
import cats.data.NonEmptyList
import cloud.nio.impl.drs.DrsCloudNioFileProvider.DrsReadInterpreter
import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, GoogleOauthDrsCredentials}
import com.google.cloud.NoCredentials
import com.typesafe.config.{Config, ConfigFactory}
import common.assertion.CromwellTimeoutSpec
import common.mock.MockSugar
import cromwell.backend.google.pipelines.common.PipelinesApiFileInput
import cromwell.backend.google.pipelines.common.io.{DiskType, PipelinesApiWorkingDisk}
import cromwell.core.path.DefaultPathBuilder
import cromwell.filesystems.drs.DrsPathBuilder
import org.mockito.Mockito._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.concurrent.duration.DurationInt

class PipelinesApiAsyncBackendJobExecutionActorSpec
    extends AnyFlatSpec
    with CromwellTimeoutSpec
    with Matchers
    with MockSugar {
  behavior of "PipelinesParameterConversions"

  it should "group files by bucket" in {

    def makeInput(bucket: String, name: String): PipelinesApiFileInput = {
      val mockCloudPath = mock[cromwell.core.path.Path]
      when(mockCloudPath.pathAsString) thenReturn s"gs://$bucket/$name"

      PipelinesApiFileInput(
        name = name,
        cloudPath = mockCloudPath,
        relativeHostPath = DefaultPathBuilder.build(Paths.get(s"$bucket/$name")),
        mount = null
      )
    }

    val inputs: List[PipelinesApiFileInput] = List(
      ("foo", "file1"),
      ("foo", "file2"),
      ("bar", "file1"),
      ("bar", "file2"),
      ("bar", "file3"),
      ("baz", "file1")
    ) map (makeInput _).tupled.apply

    val expected =
      Map("foo" -> (NonEmptyList.of(0, 1) map inputs.apply)) ++
        Map("bar" -> (NonEmptyList.of(2, 3, 4) map inputs.apply)) ++
        Map("baz" -> NonEmptyList.of(inputs(5)))

    PipelinesApiAsyncBackendJobExecutionActor.groupParametersByGcsBucket(inputs) shouldEqual expected
  }

  it should "generate a CSV manifest for DRS inputs, ignoring non-DRS inputs" in {
    def makeDrsPathBuilder: DrsPathBuilder = {
      val drsResolverConfig: Config = ConfigFactory.parseString(
        """resolver {
          |   url = "http://drshub-url"
          |}
          |""".stripMargin
      )

      val fakeCredentials = NoCredentials.getInstance

      val drsReadInterpreter: DrsReadInterpreter = (_, _) =>
        throw new UnsupportedOperationException(
          "PipelinesApiAsyncBackendJobExecutionActorSpec doesn't need to use drs read interpreter."
        )

      DrsPathBuilder(
        new DrsCloudNioFileSystemProvider(drsResolverConfig,
                                          GoogleOauthDrsCredentials(fakeCredentials, 1.minutes),
                                          drsReadInterpreter
        ),
        None
      )
    }

    val mount = PipelinesApiWorkingDisk(DiskType.LOCAL, 1)

    def makeDrsInput(name: String, drsUri: String, containerPath: String): PipelinesApiFileInput = {
      val drsPath = makeDrsPathBuilder.build(drsUri).get
      val containerRelativePath = DefaultPathBuilder.get(containerPath)
      PipelinesApiFileInput(name, drsPath, containerRelativePath, mount)
    }

    val nonDrsInput: PipelinesApiFileInput = PipelinesApiFileInput("nnn",
                                                                   DefaultPathBuilder.get("/local/nnn.bai"),
                                                                   DefaultPathBuilder.get("/path/to/nnn.bai"),
                                                                   mount
    )

    val inputs = List(
      makeDrsInput("aaa", "drs://drs.example.org/aaa", "path/to/aaa.bai"),
      nonDrsInput,
      makeDrsInput("bbb", "drs://drs.example.org/bbb", "path/to/bbb.bai")
    )

    PipelinesApiAsyncBackendJobExecutionActor.generateDrsLocalizerManifest(inputs) shouldEqual
      "drs://drs.example.org/aaa,/cromwell_root/path/to/aaa.bai\r\ndrs://drs.example.org/bbb,/cromwell_root/path/to/bbb.bai\r\n"
  }
}
