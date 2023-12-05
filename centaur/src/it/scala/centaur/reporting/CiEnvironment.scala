package centaur.reporting

import scala.util.Try

/**
  * Scala representation of the CI environment values defined in `test.inc.sh`.
  */
case class CiEnvironment(
  isCi: Option[Boolean],
  `type`: Option[String],
  branch: Option[String],
  event: Option[String],
  tag: Option[String],
  number: Option[String],
  provider: Option[String],
  url: Option[String],
  os: Option[String],
  centaurType: Option[String]
)

object CiEnvironment {
  def apply(): CiEnvironment =
    new CiEnvironment(
      isCi = sys.env.get("CROMWELL_BUILD_IS_CI").flatMap(tryToBoolean),
      `type` = sys.env.get("CROMWELL_BUILD_TYPE"),
      branch = sys.env.get("CROMWELL_BUILD_BRANCH"),
      event = sys.env.get("CROMWELL_BUILD_EVENT"),
      tag = sys.env.get("CROMWELL_BUILD_TAG"),
      number = sys.env.get("CROMWELL_BUILD_NUMBER"),
      provider = sys.env.get("CROMWELL_BUILD_PROVIDER"),
      os = sys.env.get("CROMWELL_BUILD_OS"),
      url = sys.env.get("CROMWELL_BUILD_URL"),
      centaurType = sys.env.get("CROMWELL_BUILD_CENTAUR_TYPE")
    )

  /** Try converting the value to a boolean, or return None. */
  private def tryToBoolean(string: String): Option[Boolean] = Try(string.toBoolean).toOption
}
