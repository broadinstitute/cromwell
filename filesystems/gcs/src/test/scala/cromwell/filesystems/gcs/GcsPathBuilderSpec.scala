package cromwell.filesystems.gcs

import com.google.api.gax.retrying.RetrySettings
import com.google.auth.Credentials
import com.google.auth.oauth2.ServiceAccountCredentials
import com.google.cloud.NoCredentials
import com.google.cloud.storage.contrib.nio.{CloudStorageConfiguration, CloudStorageFileSystemProvider}
import cromwell.cloudsupport.gcp.auth.ServiceAccountTestSupport
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.core.path._
import cromwell.core.{TestKitSuite, WorkflowOptions}
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.Tables.Table

import java.io.ByteArrayInputStream

class GcsPathBuilderSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with PathBuilderSpecUtils
    with ServiceAccountTestSupport {

  behavior of "GcsPathBuilder"

  it should "use google project credentials when provided in the workflow options" in {
    val wfOptionsWithProject = WorkflowOptions.fromMap(Map("google_project" -> "my_project")).get

    val gcsPathBuilderWithProjectInfo = MockGcsPathBuilder.withOptions(wfOptionsWithProject)

    gcsPathBuilderWithProjectInfo.projectId shouldBe "my_project"
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
      getNameCount = 3,
      isAbsolute = true
    ),
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
      getNameCount = 3,
      isAbsolute = true
    ),
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
      getNameCount = 3,
      isAbsolute = true
    ),
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
      getNameCount = 1,
      isAbsolute = false
    ),
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
      getNameCount = 0,
      isAbsolute = true
    ),
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
      getNameCount = 1,
      isAbsolute = true
    ),
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
      getNameCount = 2,
      isAbsolute = true
    ),

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
      getNameCount = 1,
      isAbsolute = true
    ),
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
      getNameCount = 1,
      isAbsolute = true
    ),
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
      getNameCount = 3,
      isAbsolute = true
    ),
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
      getNameCount = 3,
      isAbsolute = true
    ),

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
      getNameCount = 0,
      isAbsolute = true
    ),
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
      getNameCount = 1,
      isAbsolute = false
    ),
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
      getNameCount = 2,
      isAbsolute = true
    ),
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
      getNameCount = 1,
      isAbsolute = true
    ),
    GoodPath(
      description = "a bucket with an underscore",
      path = s"gs://hello_underscore/world",
      normalize = true,
      pathAsString = s"gs://hello_underscore/world",
      pathWithoutScheme = s"hello_underscore/world",
      parent = s"gs://hello_underscore/",
      getParent = s"gs://hello_underscore/",
      root = s"gs://hello_underscore/",
      name = "world",
      getFileName = s"gs://hello_underscore/world",
      getNameCount = 1,
      isAbsolute = true
    )
  )

  private def badPaths = Seq(
    BadPath("an empty path", "", "Path \"\" does not have a gcs scheme"),
    BadPath("an bucketless path",
            "gs://",
            "The specified GCS path 'gs://' does not parse as a URI.\nExpected authority at index 5: gs://"
    ),
    BadPath(
      "a bucket named .",
      "gs://./hello/world",
      "The path 'gs://./hello/world' does not seem to be a valid GCS path. Please check that it starts with gs:// and that the bucket and object follow GCS naming guidelines at https://cloud.google.com/storage/docs/naming."
    ),
    BadPath(
      "a non ascii bucket name",
      "gs://nonasciibucket£€/hello/world",
      "The path 'gs://nonasciibucket£€/hello/world' does not seem to be a valid GCS path. Please check that it starts with gs:// and that the bucket and object follow GCS naming guidelines at https://cloud.google.com/storage/docs/naming."
    ),
    BadPath("a https path", "https://hello/world", "Cloud Storage URIs must have 'gs' scheme: https://hello/world"),
    BadPath("a file uri path", "file:///hello/world", "Cloud Storage URIs must have 'gs' scheme: file:///hello/world"),
    BadPath("a relative file path", "hello/world", "Path \"hello/world\" does not have a gcs scheme"),
    BadPath("an absolute file path", "/hello/world", "Path \"/hello/world\" does not have a gcs scheme")
  )

  private lazy val pathBuilder =
    MockGcsPathBuilder.instance

  it should "not mix up credentials" in {
    def retrySettings: RetrySettings = RetrySettings.newBuilder().build()
    val cloudStorageConfig: CloudStorageConfiguration = GcsStorage.DefaultCloudStorageConfiguration

    val noCredentials = NoCredentials.getInstance()
    val noCredentialsPathBuilder: GcsPathBuilder = {
      val noCredentialsStorage = GcsStorage.gcsStorage("no-credentials", noCredentials, retrySettings)
      val noCredentialsStorageOptions =
        GcsStorage.gcsStorageOptions(noCredentials, retrySettings, Option("proj-no-credentials"))

      new GcsPathBuilder(noCredentialsStorage, cloudStorageConfig, noCredentialsStorageOptions)
    }

    val serviceAccountCredentials =
      ServiceAccountCredentials.fromStream(new ByteArrayInputStream(serviceAccountJsonContents.getBytes))
    val serviceAccountPathBuilder: GcsPathBuilder = {
      val serviceAccountStorage = GcsStorage.gcsStorage("service-account", serviceAccountCredentials, retrySettings)
      val serviceAccountStorageOptions =
        GcsStorage.gcsStorageOptions(serviceAccountCredentials, retrySettings, Option("proj-service-account"))

      new GcsPathBuilder(serviceAccountStorage, cloudStorageConfig, serviceAccountStorageOptions)
    }

    def credentialsForPath(gcsPath: GcsPath): Credentials = {
      val cloudFilesystemProvider =
        gcsPath.nioPath.getFileSystem.provider().asInstanceOf[CloudStorageFileSystemProvider]
      val storageOptionsField = cloudFilesystemProvider.getClass.getDeclaredField("storageOptions")
      storageOptionsField.setAccessible(true)
      val storageOptions = storageOptionsField.get(cloudFilesystemProvider)
      val credentialsField = storageOptions.getClass.getSuperclass.getSuperclass.getDeclaredField("credentials")
      credentialsField.setAccessible(true)
      credentialsField.get(storageOptions).asInstanceOf[Credentials]
    }

    val noCredentialsPath = noCredentialsPathBuilder.build("gs://no-credentials").get
    val serviceAccountPath = serviceAccountPathBuilder.build("gs://service-account").get

    credentialsForPath(noCredentialsPath) shouldBe noCredentials
    credentialsForPath(serviceAccountPath) shouldBe serviceAccountCredentials
  }
}
