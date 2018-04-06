package cromwell.engine.language

import java.util.Map.Entry

import com.typesafe.config.ConfigFactory
import cromwell.engine.language.CromwellLanguages.{CromwellLanguageName, CromwellLanguageVersion}

import scala.collection.JavaConverters._

final case class LanguageConfigurationEntry(name: CromwellLanguageName, versions: Map[CromwellLanguageVersion, LanguageConfigurationEntryFields])
final case class LanguageConfigurationEntryFields(className: String, config: Map[String, Any])

object LanguageConfiguration {
  private val LanguagesConfig = ConfigFactory.load.getConfig("languages")
  private val LanguageNames: Set[String] = LanguagesConfig.entrySet().asScala.map(findFirstKey).toSet

  val AllLanguageEntries: List[LanguageConfigurationEntry] = LanguageNames.toList map { languageName =>

    val languageConfig = LanguagesConfig.getConfig(languageName)
    val versionSet = languageConfig.getConfig("versions")
    val languageVersionNames: Set[String] = versionSet.entrySet().asScala.map(findFirstKey).toSet

    val versions = (languageVersionNames.toList map { languageVersionName =>
      val configEntry = versionSet.getConfig(s""""$languageVersionName"""")
      val className: String = configEntry.getString("language-factory")
      val factoryConfig: Map[String, Any] = if (configEntry.hasPath("config")) { configEntry.getObject("config").unwrapped().asScala.toMap } else Map.empty[String, Any]
      val fields = LanguageConfigurationEntryFields(className, factoryConfig)
      languageVersionName -> fields
    }).toMap

    LanguageConfigurationEntry(languageName, versions)
  }

  // Gets the first key in a hocon key entry (which might contain several keys, perhaps surrounded by quotes)
  // EG
  // foo {
  //   "1.0" {
  //     "x.y".z = "hello"
  //   }
  // }
  // Called on the entry on line 1, returns 'foo'
  // Called on the entry on line 2, returns '1.0'
  // Called on the entry on line 3, returns 'x.y'
  private def findFirstKey(entry: Entry[String, _]): String = {
    val NoQuoteFirstKey = """([^.]*)\.(.*)""".r
    val SingleQuotedFirstKey = """"(.*)"\.(.*)""".r

    entry.getKey match {
      case SingleQuotedFirstKey(firstKey, _) => firstKey
      case NoQuoteFirstKey(firstKey, _) => firstKey
    }
  }
}
