package cromwell.backend.impl

import java.nio.file.{FileSystem, Path}

import com.typesafe.config.Config
import cromwell.backend.async.{ExecutionHandle, NonRetryableExecution}
import cromwell.backend.impl.jes.Run.RunStatus
import cromwell.backend.{BackendJobDescriptor, BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.PathFactory._
import cromwell.core.{CallContext, WorkflowOptions}
import cromwell.filesystems.gcs.GoogleAuthMode.GoogleAuthOptions
import cromwell.filesystems.gcs.{GcsFileSystem, GcsFileSystemProvider, GoogleAuthMode, GoogleConfiguration}

import scala.util.Try


package object jes {

  implicit class PathString(val str: String) extends AnyVal {
    def isGcsUrl: Boolean = str.startsWith("gs://")
    def isUriWithProtocol: Boolean = "^[a-z]+://".r.findFirstIn(str).nonEmpty

    def toPath(fss: List[FileSystem]): Path = buildPath(str, fss)
    def toPath(fs: FileSystem): Path = str.toPath(List(fs))

    def toAbsolutePath(fss: List[FileSystem]): Path = str.toPath(fss).toAbsolutePath
    def toAbsolutePath(fs: FileSystem): Path = str.toAbsolutePath(List(fs))

    def toDirectory(fss: List[FileSystem]): Path = buildPathAsDirectory(str, fss)
    def toDirectory(fs: FileSystem): Path = str.toDirectory(List(fs))

    // TODO this needs to go away because it's gcs specific. Replacing gcs FS with google implementatio (when available) will take care of it
    private def buildPathAsDirectory(rawString: String, fileSystems: List[FileSystem]): Path = {
      findFileSystem(rawString, fileSystems, {
        case fs: GcsFileSystem => Try(fs.getPathAsDirectory(rawString))
        case fs => Try(fs.getPath(rawString))
      })
    }
  }

  implicit class EnhancedPath(val path: Path) extends AnyVal {
    def asDirectory = path.toString.toDirectory(path.getFileSystem)
  }

  /**
    * Representing a running JES execution, instances of this class are never Done and it is never okay to
    * ask them for results.
    */
  case class JesPendingExecutionHandle(jobDescriptor: BackendJobDescriptor,
                                       jesOutputs: Seq[JesFileOutput],
                                       run: Run,
                                       previousStatus: Option[RunStatus]) extends ExecutionHandle {
    override val isDone = false
    override val result = NonRetryableExecution(new IllegalStateException("JesPendingExecutionHandle cannot yield a result"))
  }

  implicit class GoogleAuthWorkflowOptions(val workflowOptions: WorkflowOptions) extends AnyVal {
    def toGoogleAuthOptions: GoogleAuthMode.GoogleAuthOptions = new GoogleAuthOptions {
      override def get(key: String): Try[String] = workflowOptions.get(key)
    }
  }

  def buildGcsFileSystem(backendConfig: Config, workflowDescriptor: BackendWorkflowDescriptor): GcsFileSystem = {
      // PBE what follows lacks any semblance of error checking

      // PBE check the config sanity in the actory factory, preferably in a fail-fast manner.
      val genomicsAuth: GoogleAuthMode =
        GoogleConfiguration.Instance.auth(backendConfig.getString("genomics.auth")) getOrElse { throw new RuntimeException("borked config") }

      // PBE It might be nice to only build this filesystems once per workflow (maybe in the initialization actor) and
      // somehow make that available to the workflow's job execution actors.
      val authOptions = workflowDescriptor.workflowOptions.toGoogleAuthOptions
      GcsFileSystemProvider(genomicsAuth.buildStorage(authOptions)).getFileSystem
  }

  object JesCallPaths {
    def apply(jobKey: BackendJobDescriptorKey, gcsFileSystem: GcsFileSystem, workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config): JesCallPaths = {
      val GcsRootOptionKey = "jes_gcs_root"
      val CallPrefix = "call"
      val ShardPrefix = "shard"
      val AttemptPrefix = "attempt"

      def jesLogBasename = {
        val index = jobKey.index.map(s => s"-$s").getOrElse("")
        s"${jobKey.scope.unqualifiedName}$index"
      }

      val rootPath: Path =
        gcsFileSystem.getPath(workflowDescriptor.workflowOptions.getOrElse(GcsRootOptionKey, backendConfig.getString("root")))

      val workflowRootPath: Path = rootPath.resolve(workflowDescriptor.workflowNamespace.workflow.unqualifiedName)
            .resolve(workflowDescriptor.id.toString)

      val callRootPath: Path = {
        val callName = jobKey.call.fullyQualifiedName.split('.').last
        val call = s"$CallPrefix-$callName"
        val shard = jobKey.index map { s => s"$ShardPrefix-$s" } getOrElse ""
        val retry = if (jobKey.attempt > 1) s"$AttemptPrefix-${jobKey.attempt}" else ""

        List(call, shard, retry).foldLeft(workflowRootPath)((path, dir) => path.resolve(dir))
      }

      val returnCodeFilename: String = s"$jesLogBasename-rc.txt"
      val stdoutFilename: String = s"$jesLogBasename-stdout.log"
      val stderrFilename: String = s"$jesLogBasename-stderr.log"
      val jesLogFilename: String = s"$jesLogBasename.log"
      JesCallPaths(rootPath, workflowRootPath, callRootPath, stdoutFilename, stderrFilename, jesLogFilename, returnCodeFilename)
    }
  }

  case class JesCallPaths private(rootPath: Path, workflowRootPath: Path, callRootPath: Path, stdoutFilename: String, stderrFilename: String, jesLogFilename: String, returnCodeFilename: String) {
    lazy val returnCodePath: Path = callRootPath.resolve(returnCodeFilename)
    lazy val stdoutPath: Path = callRootPath.resolve(stdoutFilename)
    lazy val stderrPath: Path = callRootPath.resolve(stderrFilename)

    def buildCallContext: CallContext = {
      new CallContext(callRootPath.toString, stdoutFilename, stderrFilename)
    }
  }
}
