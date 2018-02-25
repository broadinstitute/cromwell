package cromwell.backend.impl.bcs

import scala.util.{Failure, Success, Try}

trait BcsDocker {
  val image: String
}

final case class BcsDockerWithoutPath(image: String) extends BcsDocker
final case class BcsDockerWithPath(image: String, path: String) extends BcsDocker


object BcsDocker{
  val dockerWithPathPattern = s"""(\\S+)\\s+(\\S+)""".r
  val dockerWithoutPathPatter = s"""(\\S+)""".r

  def parse(s: String): Try[BcsDocker] = {
    s match {
      case dockerWithoutPathPatter(dockerImage) => Success(BcsDockerWithoutPath(dockerImage))
      case dockerWithPathPattern(dockerImage, dockerPath) => Success(BcsDockerWithPath(dockerImage, dockerPath))
      case _ => Failure(new IllegalArgumentException("must be 'ubuntu/latest oss://docker-reg/ubuntu/'"))
    }
  }
}
