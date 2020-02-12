package cromwell.backend.google.pipelines.v2beta

import java.nio.file.Paths

import cats.data.NonEmptyList
import cromwell.backend.google.pipelines.common.PipelinesApiFileInput
import cromwell.core.path.DefaultPathBuilder
import org.mockito.Mockito._
import org.scalatest.mockito.MockitoSugar
import org.scalatest.{FlatSpec, Matchers}

class PipelinesApiAsyncBackendJobExecutionActorSpec extends FlatSpec with Matchers with MockitoSugar {
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
}
