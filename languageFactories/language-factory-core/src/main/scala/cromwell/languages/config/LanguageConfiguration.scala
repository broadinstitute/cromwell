package cromwell.languages.config

import java.util.Map.Entry

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.languages.config.CromwellLanguages.{CromwellLanguageName, CromwellLanguageVersion}

import scala.collection.JavaConverters._

final case class LanguagesConfiguration(languages: List[LanguageVersionConfigurationEntry], default: Option[String])
final case class LanguageVersionConfigurationEntry(name: CromwellLanguageName, versions: Map[CromwellLanguageVersion, LanguageVersionConfig], default: Option[String])
final case class LanguageVersionConfig(className: String, config: Config)

object LanguageConfiguration {
  private val LanguagesConfig = ConfigFactory.load.getConfig("languages")
  private val DefaultLanguageName: Option[String] = if (LanguagesConfig.hasPath("default")) Option(LanguagesConfig.getString("default")) else None

  private val LanguageNames: Set[String] = LanguagesConfig.entrySet().asScala.map(findFirstKey).filterNot(_ == "default").toSet

  val AllLanguageEntries: LanguagesConfiguration = {
    val languages = LanguageNames.toList map { languageName =>

      val languageConfig = LanguagesConfig.getConfig(languageName)
      val defaultVersionName: Option[String] = if (LanguagesConfig.hasPath("default")) { Option(LanguagesConfig.getString("default")) } else None
      val versionSet = languageConfig.getConfig("versions")
      val languageVersionNames: Set[String] = versionSet.entrySet().asScala.map(findFirstKey).filterNot(_ == "default").toSet

      val versions = (languageVersionNames.toList map { languageVersionName =>
        val configEntry = versionSet.getConfig(s""""$languageVersionName"""")
        val className: String = configEntry.getString("language-factory")
        val factoryConfig: Config = if (configEntry.hasPath("config")) configEntry.getConfig("config") else ConfigFactory.empty()
        val fields = LanguageVersionConfig(className, factoryConfig)
        languageVersionName -> fields
      }).toMap

      LanguageVersionConfigurationEntry(languageName, versions, defaultVersionName)
    }

    LanguagesConfiguration(languages, DefaultLanguageName)
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
    val NoQuoteFirstKey = """([^/.]*)\.(.*)""".r
    val SingleQuotedFirstKey = """"(.*)"\.(.*)""".r
    val NoQuoteNoDot = """([^\.\"]*)""".r

    entry.getKey match {
      case SingleQuotedFirstKey(firstKey, _) => firstKey
      case NoQuoteFirstKey(firstKey, _) => firstKey
      case NoQuoteNoDot(key) => key
      case other => throw new IllegalArgumentException(s"Unexpected key format: $other")
    }
  }
}
