package cromwell.engine.language

import cromwell.engine.language.CromwellLanguages._
import cromwell.languages.LanguageFactory

// Construct a singleton instance of this class using 'initLanguages' below.
final case class CromwellLanguages private(languageConfig: List[LanguageConfigurationEntry]) {

  val languages: Map[CromwellLanguageName, LanguageVersions] = makeLanguages

  private def makeLanguages = (languageConfig map { lc =>
    val versions = lc.versions map { case (languageVersion, languageConfigEntryFields) =>
      languageVersion -> makeLanguageFactory(languageConfigEntryFields.className, languageConfigEntryFields.config)
    }

    lc.name.toUpperCase -> LanguageVersions(versions)
  }).toMap


  private def makeLanguageFactory(className: String, config: Map[String, Any]) = {
    Class.forName(className)
      .getConstructor(classOf[Map[String, Any]])
      .newInstance(config)
      .asInstanceOf[LanguageFactory]
  }
}

/**
  * Holds all the registered versions of a language.
  */
final case class LanguageVersions private(allVersions: Map[CromwellLanguageVersion, LanguageFactory])

object CromwellLanguages {
  type CromwellLanguageName = String
  type CromwellLanguageVersion = String

  private var _instance: CromwellLanguages = _
  lazy val instance: CromwellLanguages = _instance

  def initLanguages(backendEntries: List[LanguageConfigurationEntry]): Unit = {
    _instance = CromwellLanguages(backendEntries)
  }
}
