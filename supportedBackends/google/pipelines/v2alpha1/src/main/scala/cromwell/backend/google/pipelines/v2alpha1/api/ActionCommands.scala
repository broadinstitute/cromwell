package cromwell.backend.google.pipelines.v2alpha1.api
import akka.http.scaladsl.model.ContentType
import common.util.StringUtil._
import cromwell.core.path.Path
import mouse.all._
import org.apache.commons.text.StringEscapeUtils.ESCAPE_XSI

import scala.concurrent.duration._

/**
  * Utility methods to build shell commands for localization / delocalization.
  */
object ActionCommands {
  implicit val nbRetries: Int = 5
  implicit val waitBetweenRetries: FiniteDuration = 5.seconds

  implicit class ShellPath(val path: Path) extends AnyVal {
    // The command String runs in Bourne shell so shell metacharacters in filenames must be escaped
    def escape = ESCAPE_XSI.translate(path.pathAsString)
  }

  def makeContainerDirectory(containerPath: Path) = s"mkdir -p ${containerPath.escape}"

  def delocalizeDirectory(containerPath: Path, cloudPath: Path) = retry {
    s"gsutil -m rsync -r ${containerPath.escape} ${cloudPath.escape}"
  }

  /**
    * As per https://cloud.google.com/storage/docs/gsutil/addlhelp/HowSubdirectoriesWork, rule #2
    * If one attempts a 
    *  gsutil cp /local/file.txt gs://bucket/subdir/file.txt
    *  AND 
    *  there exists a folder gs://bucket/subdir/file.txt_thisCouldBeAnything
    *  then gs://bucket/subdir/file.txt will be treated as a directory, and /local/file.txt will be copied under gs://bucket/subdir/file.txt/file.txt
    *  and not gs://bucket/subdir/file.txt.
    *
    * By instead using the parent directory (and ensuring it ends with a slash), gsutil will treat that as a directory and put the file under it.
    * So the final gsutil command will look something like gsutil cp /local/file.txt gs://bucket/subdir/
    */
  def delocalizeFile(containerPath: Path, cloudPath: Path, contentType: Option[ContentType]) = retry {
    val contentTypeFlag = contentType.map(ct => s"""-h "Content-Type: $ct"""").getOrElse("")
    s"gsutil $contentTypeFlag cp ${containerPath.escape} ${cloudPath.parent.escape.ensureSlashed}"
  }

  /**
    * delocalizeFile necessarily copies the file to the same name. Use this if you want to to specify a name different from the original
    * Make sure that there's no object named "yourfinalname_something" (see above) in the same cloud directory.
    */
  def delocalizeFileTo(containerPath: Path, cloudPath: Path, contentType: Option[ContentType]) = retry {
    val contentTypeFlag = contentType.map(ct => s"""-h "Content-Type: $ct"""").getOrElse("")
    s"gsutil $contentTypeFlag cp ${containerPath.escape} ${cloudPath.escape}"
  }

  def ifExist(containerPath: Path)(f: => String) = s"if [[ -e ${containerPath.escape} ]]; then $f; fi"

  def every(duration: FiniteDuration)(f: => String) = s"while true; do $f 2> /dev/null || true; sleep ${duration.toSeconds}; done"

  def retry(f: => String)(implicit times: Int, wait: FiniteDuration) = {
    s"for i in `seq $times`; do $f && break; sleep ${wait.toSeconds}; done"
  }

  def delocalizeFileOrDirectory(containerPath: Path, cloudPath: Path, contentType: Option[ContentType]) = {
    s"if [[ -d ${containerPath.escape} ]]; then ${delocalizeDirectory(containerPath, cloudPath)}; else ${delocalizeFile(containerPath, cloudPath, contentType)}; fi"
  }

  def localizeDirectory(cloudPath: Path, containerPath: Path) = retry {
    s"${containerPath |> makeContainerDirectory} && gsutil -m rsync -r ${cloudPath.escape} ${containerPath.escape}"
  }

  def localizeFile(cloudPath: Path, containerPath: Path) = retry {
    s"gsutil cp ${cloudPath.escape} ${containerPath.escape}"
  }
}
