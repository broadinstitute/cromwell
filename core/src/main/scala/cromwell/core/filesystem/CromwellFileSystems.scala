package cromwell.core.filesystem

import java.lang.reflect.Constructor

import cats.instances.list._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.typesafe.config.{Config, ConfigFactory, ConfigObject}
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.core.path.{DefaultPathBuilderFactory, PathBuilderFactory}
import net.ceedubs.ficus.Ficus._
import shapeless.syntax.typeable._

import scala.collection.JavaConverters._
import scala.util.{Failure, Try}

/**
  * Validates filesystem configuration and provide methods to build PathBuilderFactories.
  */
class CromwellFileSystems(globalConfig: Config) {
  // Validate the configuration and creates a Map of PathBuilderFactory constructors
  private [filesystem] val factoryBuilders: Map[String, Constructor[_]] = if (globalConfig.hasPath("filesystems")) {
    val rawConfigSet = globalConfig.getObject("filesystems").entrySet.asScala
    val configMap = rawConfigSet.toList.map({ entry => entry.getKey -> entry.getValue })
    val constructorMap = configMap.traverse[ErrorOr, (String, Constructor[_])]({
      case (key, fsConfig: ConfigObject) => createConstructor(key, fsConfig).toValidated
      case (key, _) => s"Invalid filesystem configuration for $key".invalidNel
    }).map(_.toMap)

    constructorMap.unsafe("Failed to initialize Cromwell filesystems")
  } else Map.empty

  // Create a constructor from a configuration object
  private def createConstructor(key: String, configObject: ConfigObject): Checked[(String, Constructor[_])] = for {
    clazz <- configObject.toConfig.as[Option[String]]("class").toChecked(s"Filesystem configuration $key doesn't have a class field")
    constructor <- createConstructor(key, clazz)
  } yield constructor

  // Create a constructor from a class name
  private def createConstructor(filesystem: String, className: String): Checked[(String, Constructor[_])] = Try (
    filesystem -> Class.forName(className).getConstructor(classOf[Config], classOf[Config])
  ).recoverWith({
    case e: ClassNotFoundException => Failure(
      new RuntimeException(s"Class $className for filesystem $filesystem cannot be found in the class path.", e)
    )
    case e: NoSuchMethodException => Failure(
      new RuntimeException(s"Class $className for filesystem $filesystem does not have the required constructor signature.", e)
    )
  }).toChecked

  // Instantiate a PathBuilderFactor from its constructor and instance config
  private def instantiate(name: String, constructor: Constructor[_], instanceConfig: Config): Checked[PathBuilderFactory] = {
    for {
      instance <- Try(constructor.newInstance(globalConfig, instanceConfig)).toChecked
      cast <- instance.cast[PathBuilderFactory].toChecked(s"The filesystem class for $name is not an instance of PathBuilderFactory")
    } yield cast
  }

  // Look for a constructor in the map of known filesystems
  private def getConstructor(fileSystemName: String): Checked[Constructor[_]] = factoryBuilders
    .get(fileSystemName)
    .toChecked(s"Cannot find a filesystem with name $fileSystemName in the configuration. Available filesystems: ${factoryBuilders.keySet.mkString(", ")}")

  /**
    * Try to find a configured filesystem with the given name and to build a PathFactory for it
    * @param name name of the filesystem
    * @param instanceConfig filesystem specific configuration for this instance of the factory to build
    */
  def buildFactory(name: String, instanceConfig: Config): Checked[PathBuilderFactory] = {
    if (DefaultPathBuilderFactory.name.equalsIgnoreCase(name)) DefaultPathBuilderFactory.validNelCheck
    else for {
      constructor <- getConstructor(name)
      factory <- instantiate(name, constructor, instanceConfig)
    } yield factory
  }

  /**
    * Given a filesystems config, build the PathBuilderFactories
    */
  def factoriesFromConfig(filesystemsConfig: Config): Checked[Map[String, PathBuilderFactory]] = {
    if (filesystemsConfig.hasPath("filesystems")) {
      // Iterate over the config entries under the "filesystems" config
      val rawConfigSet = filesystemsConfig.getObject("filesystems").entrySet().asScala
      val configMap = rawConfigSet.toList.map({ entry => entry.getKey -> entry.getValue })

      configMap.traverse[ErrorOr, (String, PathBuilderFactory)]({
        // build the factory for each entry
        case (key, config: ConfigObject) => buildFactory(key, config.toConfig).toValidated.map(key -> _)
        case (key, _) => s"Invalid filesystem backend configuration for $key".invalidNel
      }).map(_.toMap).toEither
    } else Map.empty[String, PathBuilderFactory].validNelCheck
  }
}

object CromwellFileSystems {
  val instance = new CromwellFileSystems(ConfigFactory.load())
}
