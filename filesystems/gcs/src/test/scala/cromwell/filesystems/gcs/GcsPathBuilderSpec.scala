package cromwell.filesystems.gcs

import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration
import cromwell.core.path._
import cromwell.core.{TestKitSuite, WorkflowOptions}
import cromwell.filesystems.gcs.auth.{GoogleAuthMode, GoogleAuthModeSpec}
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpecLike, Matchers}

class GcsPathBuilderSpec extends TestKitSuite with FlatSpecLike with Matchers with PathBuilderSpecUtils {

  behavior of "GcsPathBuilder"

  it should "use google project credentials when provided in the workflow options" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()

    val wfOptionsWithProject = WorkflowOptions.fromMap(Map("google_project" -> "my_project")).get

    val gcsPathBuilderWithProjectInfo = new GcsPathBuilder(
      GoogleAuthMode.MockAuthMode,
      "cromwell-test",
      None,
      CloudStorageConfiguration.DEFAULT,
      wfOptionsWithProject
    )

    gcsPathBuilderWithProjectInfo.getProjectId shouldBe "my_project"
  }

  it should behave like truncateCommonRoots(pathBuilder, pathsToTruncate)

  goodPaths foreach { goodPath =>
    it should behave like buildGoodPath(pathBuilder, goodPath)
  }

  badPaths foreach { badPath =>
    it should behave like buildBadPath(pathBuilder, badPath)
  }

  private def pathsToTruncate = Table(
    ("context", "file", "relative"),
    ("gs://bucket", "gs://bucket/path/to/file", "path/to/file"),
    ("gs://bucket/path/to/my/dir", "gs://bucket/path/to/my/dir/file", "file"),
    ("gs://bucket/path/to/my/dir", "gs://bucket/path/to/my/dir//file", "file"),
    // NOTE: Next two are different from the DefaultPathBuilder. "//" doesn't build to "/" in the GcsPathBuilder
    ("gs://bucket/path/to/my//dir", "gs://bucket/path/to/my/dir/file", "dir/file"),
    ("gs://bucket/path/to/my//dir", "gs://bucket/path/to/my/dir//file", "dir//file"),
    ("gs://bucket/path/to/my/dir", "gs://bucket/path/./to/my/dir/file", "./to/my/dir/file"),
    ("gs://bucket/path/to/my/dir/with/file", "gs://bucket/path/to/other/dir/with/file", "other/dir/with/file")
  )

  private def bucket = "mymadeupbucket"

  private def goodPaths = Seq(
    GoodPath(
      description = "a path with spaces",
      path = s"gs://$bucket/hello/world/with spaces",
      normalize = false,
      pathAsString = s"gs://$bucket/hello/world/with spaces",
      pathWithoutScheme = s"$bucket/hello/world/with spaces",
      parent = s"gs://$bucket/hello/world/",
      getParent = s"gs://$bucket/hello/world/",
      root = s"gs://$bucket/",
      name = "with spaces",
      getFileName = s"gs://$bucket/with spaces",
      toUriHost = bucket,
      toUriPath = "/hello/world/with%20spaces",
      toUriStartsWith = s"gs://$bucket/hello/world/with%20spaces",
      toUriEndsWith = s"gs://$bucket/hello/world/with%20spaces",
      getNameCount = 3,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "a path with non-ascii",
      path = s"gs://$bucket/hello/world/with non ascii £€",
      normalize = false,
      pathAsString = s"gs://$bucket/hello/world/with non ascii £€",
      pathWithoutScheme = s"$bucket/hello/world/with non ascii £€",
      parent = s"gs://$bucket/hello/world/",
      getParent = s"gs://$bucket/hello/world/",
      root = s"gs://$bucket/",
      name = "with non ascii £€",
      getFileName = s"gs://$bucket/with non ascii £€",
      toUriHost = bucket,
      toUriPath = "/hello/world/with%20non%20ascii%20£€",
      toUriStartsWith = s"gs://$bucket/hello/world/with%20non%20ascii%20£€",
      toUriEndsWith = s"gs://$bucket/hello/world/with%20non%20ascii%20£€",
      getNameCount = 3,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "a gs uri path with encoded characters",
      path = s"gs://$bucket/hello/world/encoded%20spaces",
      normalize = false,
      pathAsString = s"gs://$bucket/hello/world/encoded%20spaces",
      pathWithoutScheme = s"$bucket/hello/world/encoded%20spaces",
      parent = s"gs://$bucket/hello/world/",
      getParent = s"gs://$bucket/hello/world/",
      root = s"gs://$bucket/",
      name = "encoded%20spaces",
      getFileName = s"gs://$bucket/encoded%20spaces",
      toUriHost = bucket,
      toUriPath = "/hello/world/encoded%2520spaces",
      toUriStartsWith = s"gs://$bucket/hello/world/encoded%2520spaces",
      toUriEndsWith = s"gs://$bucket/hello/world/encoded%2520spaces",
      getNameCount = 3,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "a bucket only path",
      path = s"gs://$bucket",
      normalize = false,
      pathAsString = s"gs://$bucket/",
      pathWithoutScheme = s"$bucket/",
      parent = null,
      getParent = null,
      root = s"gs://$bucket/",
      name = "",
      getFileName = s"gs://$bucket/",
      toUriHost = bucket,
      toUriPath = "/",
      toUriStartsWith = s"gs://$bucket/",
      toUriEndsWith = s"gs://$bucket/",
      getNameCount = 1,
      isAbsolute = false,
      isDirectory = true),

    GoodPath(
      description = "a bucket only path ending in a /",
      path = s"gs://$bucket/",
      normalize = false,
      pathAsString = s"gs://$bucket/",
      pathWithoutScheme = s"$bucket/",
      parent = null,
      getParent = null,
      root = s"gs://$bucket/",
      name = "",
      getFileName = null,
      toUriHost = bucket,
      toUriPath = "/",
      toUriStartsWith = s"gs://$bucket/",
      toUriEndsWith = s"gs://$bucket/",
      getNameCount = 0,
      isAbsolute = true,
      isDirectory = true),

    GoodPath(
      description = "a file at the top of the bucket",
      path = s"gs://$bucket/hello",
      normalize = false,
      pathAsString = s"gs://$bucket/hello",
      pathWithoutScheme = s"$bucket/hello",
      parent = s"gs://$bucket/",
      getParent = s"gs://$bucket/",
      root = s"gs://$bucket/",
      name = "hello",
      getFileName = s"gs://$bucket/hello",
      toUriHost = bucket,
      toUriPath = "/hello",
      toUriStartsWith = s"gs://$bucket/hello",
      toUriEndsWith = s"gs://$bucket/hello",
      getNameCount = 1,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "a path ending in /",
      path = s"gs://$bucket/hello/world/",
      normalize = false,
      pathAsString = s"gs://$bucket/hello/world/",
      pathWithoutScheme = s"$bucket/hello/world/",
      parent = s"gs://$bucket/hello/",
      getParent = s"gs://$bucket/hello/",
      root = s"gs://$bucket/",
      name = "world",
      getFileName = s"gs://$bucket/world",
      toUriHost = bucket,
      toUriPath = "/hello/world/",
      toUriStartsWith = s"gs://$bucket/hello/world/",
      toUriEndsWith = s"gs://$bucket/hello/world/",
      getNameCount = 2,
      isAbsolute = true,
      isDirectory = true),

    // Special paths

    GoodPath(
      description = "a bucket with a path .",
      path = s"gs://$bucket/.",
      normalize = false,
      pathAsString = s"gs://$bucket/.",
      pathWithoutScheme = s"$bucket/.",
      parent = null,
      getParent = s"gs://$bucket/",
      root = s"gs://$bucket/",
      name = "",
      getFileName = s"gs://$bucket/.",
      toUriHost = bucket,
      toUriPath = "/.",
      toUriStartsWith = s"gs://$bucket/.",
      toUriEndsWith = s"gs://$bucket/.",
      getNameCount = 1,
      isAbsolute = true,
      isDirectory = true),

    GoodPath(
      description = "a bucket with a path ..",
      path = s"gs://$bucket/..",
      normalize = false,
      pathAsString = s"gs://$bucket/..",
      pathWithoutScheme = s"$bucket/..",
      parent = null,
      getParent = s"gs://$bucket/",
      root = null,
      name = "",
      getFileName = s"gs://$bucket/..",
      toUriHost = bucket,
      toUriPath = "/..",
      toUriStartsWith = s"gs://$bucket/..",
      toUriEndsWith = s"gs://$bucket/..",
      getNameCount = 1,
      isAbsolute = true,
      isDirectory = true),

    GoodPath(
      description = "a bucket including . in the path",
      path = s"gs://$bucket/hello/./world",
      normalize = false,
      pathAsString = s"gs://$bucket/hello/./world",
      pathWithoutScheme = s"$bucket/hello/./world",
      parent = s"gs://$bucket/hello/",
      getParent = s"gs://$bucket/hello/./",
      root = s"gs://$bucket/",
      name = "world",
      getFileName = s"gs://$bucket/world",
      toUriHost = bucket,
      toUriPath = "/hello/./world",
      toUriStartsWith = s"gs://$bucket/hello/./world",
      toUriEndsWith = s"gs://$bucket/hello/./world",
      getNameCount = 3,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "a bucket including .. in the path",
      path = s"gs://$bucket/hello/../world",
      normalize = false,
      pathAsString = s"gs://$bucket/hello/../world",
      pathWithoutScheme = s"$bucket/hello/../world",
      parent = s"gs://$bucket/",
      getParent = s"gs://$bucket/hello/../",
      root = s"gs://$bucket/",
      name = "world",
      getFileName = s"gs://$bucket/world",
      toUriHost = bucket,
      toUriPath = "/hello/../world",
      toUriStartsWith = s"gs://$bucket/hello/../world",
      toUriEndsWith = s"gs://$bucket/hello/../world",
      getNameCount = 3,
      isAbsolute = true,
      isDirectory = false),

    // Normalized

    GoodPath(
      description = "a bucket with a normalized path .",
      path = s"gs://$bucket/.",
      normalize = true,
      pathAsString = s"gs://$bucket/",
      pathWithoutScheme = s"$bucket/",
      parent = null,
      getParent = null,
      root = s"gs://$bucket/",
      name = "",
      getFileName = null,
      toUriHost = bucket,
      toUriPath = "/",
      toUriStartsWith = s"gs://$bucket/",
      toUriEndsWith = s"gs://$bucket/",
      getNameCount = 0,
      isAbsolute = true,
      isDirectory = true),

    GoodPath(
      description = "a bucket with a normalized path ..",
      path = s"gs://$bucket/..",
      normalize = true,
      pathAsString = s"gs://$bucket/",
      pathWithoutScheme = s"$bucket/",
      parent = null,
      getParent = null,
      root = s"gs://$bucket/",
      name = "",
      getFileName = s"gs://$bucket/",
      toUriHost = bucket,
      toUriPath = "/",
      toUriStartsWith = s"gs://$bucket/",
      toUriEndsWith = s"gs://$bucket/",
      getNameCount = 1,
      isAbsolute = false,
      isDirectory = true),

    GoodPath(
      description = "a bucket including . in the normalized path",
      path = s"gs://$bucket/hello/./world",
      normalize = true,
      pathAsString = s"gs://$bucket/hello/world",
      pathWithoutScheme = s"$bucket/hello/world",
      parent = s"gs://$bucket/hello/",
      getParent = s"gs://$bucket/hello/",
      root = s"gs://$bucket/",
      name = "world",
      getFileName = s"gs://$bucket/world",
      toUriHost = bucket,
      toUriPath = "/hello/world",
      toUriStartsWith = s"gs://$bucket/hello/world",
      toUriEndsWith = s"gs://$bucket/hello/world",
      getNameCount = 2,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "a bucket including .. in the normalized path",
      path = s"gs://$bucket/hello/../world",
      normalize = true,
      pathAsString = s"gs://$bucket/world",
      pathWithoutScheme = s"$bucket/world",
      parent = s"gs://$bucket/",
      getParent = s"gs://$bucket/",
      root = s"gs://$bucket/",
      name = "world",
      getFileName = s"gs://$bucket/world",
      toUriHost = bucket,
      toUriPath = "/world",
      toUriStartsWith = s"gs://$bucket/world",
      toUriEndsWith = s"gs://$bucket/world",
      getNameCount = 1,
      isAbsolute = true,
      isDirectory = false)
  )

  private def badPaths = Seq(
    BadPath("an empty path", "", " does not have a gcs scheme"),
    BadPath("an bucketless path", "gs://", "Expected authority at index 5: gs://"),
    BadPath("a bucket named .", "gs://./hello/world", "gs://./hello/world does not have a host"),
    BadPath("a non ascii bucket name", "gs://nonasciibucket£€/hello/world",
      "gs://nonasciibucket%C2%A3%E2%82%AC/hello/world does not have a host"),
    BadPath("a https path", "https://hello/world", "Cloud Storage URIs must have 'gs' scheme: https://hello/world"),
    BadPath("a file uri path", "file:///hello/world", "Cloud Storage URIs must have 'gs' scheme: file:///hello/world"),
    BadPath("a relative file path", "hello/world", "hello/world does not have a gcs scheme"),
    BadPath("an absolute file path", "/hello/world", "/hello/world does not have a gcs scheme")
  )

  private lazy val pathBuilder = {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()

    new GcsPathBuilder(
      GoogleAuthMode.MockAuthMode,
      "cromwell-test",
      None,
      CloudStorageConfiguration.DEFAULT,
      WorkflowOptions.empty
    )
  }
}
