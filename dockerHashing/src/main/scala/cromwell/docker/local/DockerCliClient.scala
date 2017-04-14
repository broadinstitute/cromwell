package cromwell.docker.local

import scala.util.Try

/**
  * Wrapper around the docker cli.
  * https://docs.docker.com/engine/reference/commandline/docker/
  *
  * An alternative to using REST calls directly to the docker engine that has several versions, for example:
  * https://docs.docker.com/engine/api/v1.27/
  */
trait DockerCliClient {
  /**
    * Looks up a docker hash.
    *
    * @param dockerCliKey The docker hash to lookup.
    * @return The hash if found, None if not found, and Failure if an error occurs.
    */
  def lookupHash(dockerCliKey: DockerCliKey): Try[Option[String]] = {
    /*
    The stdout contains the tab separated repository/tag/digest for __all__ local images.
    Would be great to just get a single hash using the key... unfortunately
    https://github.com/docker/docker/issues/29901
     */
    forRun("docker", "images", "--digests", "--format", """{{printf "%s\t%s\t%s" .Repository .Tag .Digest}}""") {
      _.flatMap(parseHashLine).find(_.key == dockerCliKey).map(_.digest)
    }
  }

  /**
    * Pulls a docker image.
    * @param dockerCliKey The docker hash to lookup.
    * @return Failure if an error occurs.
    */
  def pull(dockerCliKey: DockerCliKey): Try[Unit] = {
    forRun("docker", "pull", dockerCliKey.fullName) { _ => () }
  }

  /**
    * Tries to run the command, then feeds the stdout to `f`. If the exit code is non-zero, returns a `Failure` with
    * a `RuntimeException` containing the stderr.
    *
    * @param cmd Command line to run.
    * @param f The function to run on the stdout contents.
    * @tparam A Return type.
    * @return An attempt to run A.
    */
  private def forRun[A](cmd: String*)(f: Seq[String] => A): Try[A] = {
    Try {
      val dockerCliResult = run(cmd)
      if (dockerCliResult.exitCode == 0) {
        f(dockerCliResult.stdout)
      } else {
        throw new RuntimeException(
          s"""|Error running: ${cmd.mkString(" ")}
              |Exit code: ${dockerCliResult.exitCode}
              |${dockerCliResult.stderr.mkString("\n")}
              |""".stripMargin)
      }
    }
  }

  /**
    * Run a command and return the result. Overridable for testing.
    *
    * @param cmd The command to run.
    * @return The results of the run wrapped in a DockerCliResult.
    */
  private[local] def run(cmd: Seq[String]): DockerCliResult = {
    import sys.process._
    var stdout = List.empty[String]
    var stderr = List.empty[String]
    val exitCode = cmd.!(ProcessLogger(line => stdout :+= line, line => stderr :+= line))
    DockerCliResult(exitCode, stdout, stderr)
  }

  /**
    * Parses a line for `lookupHash`, returning None for lines that contain the string "<none>" for any of the columns.
    * @param hashLine The line output by the stdout of the `lookupHash` command.
    * @return An optional `DockerCliHash`, if the all the columns are found.
    */
  private[local] def parseHashLine(hashLine: String): Option[DockerCliHash] = {
    val none = "<none>"
    val tokens = hashLine.split("\t").lift
    for {
      repository <- tokens(0) if repository != none
      tag <- tokens(1) if tag != none
      digest <- tokens(2) if digest != none
    } yield DockerCliHash(DockerCliKey(repository, tag), digest)
  }
}

object DockerCliClient extends DockerCliClient

/**
  * Used by `lookupHash` to return a mapping of keys to digests.
  * See note in `lookupHash` regarding why we even have to loop over all images.
  */
case class DockerCliHash(key: DockerCliKey, digest: String)

/**
  * Utility that encapsulates a command's exit code, stdout, and stderr.
  * Also used by tests to simulate cli responses.
  */
case class DockerCliResult(exitCode: Int, stdout: Seq[String], stderr: Seq[String])
