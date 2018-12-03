package cromwell.languages.config

import com.typesafe.config.Config
import cromwell.languages.config.CromwellLanguages.{CromwellLanguageName, CromwellLanguageVersion}
import cromwell.languages.LanguageFactory

// Construct a singleton instance of this class using 'initLanguages' below.
final case class CromwellLanguages private(languageConfig: LanguagesConfiguration) {

  val languages: Map[CromwellLanguageName, LanguageVersions] = makeLanguages
  val default: LanguageVersions = languages.find(lang => languageConfig.default.contains(lang._1)).getOrElse(languages.head)._2

  private def makeLanguages: Map[CromwellLanguageName, LanguageVersions] = (languageConfig.languages map { lc =>
    val versions = lc.versions map { case (languageVersion, languageConfigEntryFields) =>
      languageVersion -> makeLanguageFactory(languageConfigEntryFields.className, languageConfigEntryFields.config)
    }
    val default: LanguageFactory = versions.find(v => lc.default.contains(v._1)).getOrElse(versions.head)._2

    lc.name.toUpperCase -> LanguageVersions(versions, default)
  }).toMap

  private def makeLanguageFactory(className: String, config: Config) = {
    Class.forName(className)
      .getConstructor(classOf[Config])
      .newInstance(config)
      .asInstanceOf[LanguageFactory]
  }
}

/**
  * Holds all the registered versions of a language.
  */
final case class LanguageVersions private(allVersions: Map[CromwellLanguageVersion, LanguageFactory], default: LanguageFactory)

object CromwellLanguages {
  type CromwellLanguageName = String
  type CromwellLanguageVersion = String

  private var _instance: CromwellLanguages = _
  lazy val instance: CromwellLanguages = _instance

  def initLanguages(backendEntries: LanguagesConfiguration): Unit = {
    _instance = CromwellLanguages(backendEntries)
  }
}
