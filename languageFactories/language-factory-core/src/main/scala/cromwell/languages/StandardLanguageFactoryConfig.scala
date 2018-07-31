package cromwell.languages

final case class StandardLanguageFactoryConfig(strictValidation: Boolean, enabled: Boolean, unusedKeys: Map[String, Any])

object StandardLanguageFactoryConfig {
  def parse(values: Map[String, Any], allowExtras: Boolean): StandardLanguageFactoryConfig = {

    def parseFold(current: StandardLanguageFactoryConfig, next: (String, Any)) = next match {
      case ("strict-validation", b: Boolean) => current.copy(strictValidation = b)
      case ("enabled", b: Boolean) => current.copy(enabled = b)
      case other if allowExtras => current.copy(unusedKeys = current.unusedKeys + other)
      case other => throw new RuntimeException(s"Non-standard configuration value for language factory: $other")
    }

    val initial = StandardLanguageFactoryConfig(strictValidation = true, enabled = true, Map.empty)
    values.foldLeft(initial)(parseFold)
  }
}
