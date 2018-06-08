package cromwell.backend.google.pipelines.v2alpha1.api
import common.util.StringUtil._
import cromwell.core.path.Path
import mouse.all._
import org.apache.commons.text.StringEscapeUtils.ESCAPE_XSI

/**
  * Utility methods to build shell commands for localization / delocalization.
  */
object ActionCommands {
  implicit class ShellPath(val path: Path) extends AnyVal {
    // The command String runs in Bourne shell so shell metacharacters in filenames must be escaped
    def escape = ESCAPE_XSI.translate(path.pathAsString)
  }

  def makeContainerDirectory(containerPath: Path) = s"mkdir -p ${containerPath.escape}"

  def delocalizeDirectory(containerPath: Path, cloudPath: Path) = {
    s"gsutil -m rsync -r ${containerPath.escape} ${cloudPath.escape}"
  }

  /*
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
  def delocalizeFile(containerPath: Path, cloudPath: Path) = {
    s"gsutil cp ${containerPath.escape} ${cloudPath.parent.escape.ensureSlashed}"
  }

  def ifExist(containerPath: Path)(f: => String) = s"if [[ -e ${containerPath.escape} ]]; then $f; fi"

  def delocalizeFileOrDirectory(containerPath: Path, cloudPath: Path) = {
    s"if [[ -d ${containerPath.escape} ]]; then ${delocalizeDirectory(containerPath, cloudPath)}; else ${delocalizeFile(containerPath, cloudPath)}; fi"
  }

  def localizeDirectory(cloudPath: Path, containerPath: Path) = {
    s"${containerPath |> makeContainerDirectory} && gsutil -m rsync -r ${cloudPath.escape} ${containerPath.escape}"
  }

  def localizeFile(cloudPath: Path, containerPath: Path) = {
    s"gsutil cp ${cloudPath.escape} ${containerPath.escape}"
  }
}
