package cromwell.core.filesystem

import java.lang.reflect.Constructor

import cats.instances.list._
import cats.instances.option._
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
import scala.language.{existentials, postfixOps}
import scala.util.{Failure, Try}

/**
  * Validates filesystem configuration and provides methods to build PathBuilderFactories.
  */
class CromwellFileSystems(globalConfig: Config) {
  // Validate the configuration and creates a Map of PathBuilderFactory constructors, along with their optional singleton config
  private [filesystem] val factoryBuilders: Map[String, (Constructor[_], Option[AnyRef])] = if (globalConfig.hasPath("filesystems")) {
    val rawConfigSet = globalConfig.getObject("filesystems").entrySet.asScala
    val configMap = rawConfigSet.toList.map({ entry => entry.getKey -> entry.getValue })
    val constructorMap = configMap.traverse[ErrorOr, (String, (Constructor[_], Option[AnyRef]))]({
      case (key, fsConfig: ConfigObject) => processFileSystem(key, fsConfig)
      case (key, _) => s"Invalid filesystem configuration for $key".invalidNel
    }).map(_.toMap)

    constructorMap.unsafe("Failed to initialize Cromwell filesystems")
  } else Map.empty

  val supportedFileSystems: Iterable[String] = factoryBuilders.keys

  // Generate the appropriate constructor and optional singleton instance for a filesystem
  private def processFileSystem(key: String, fsConfig: ConfigObject): ErrorOr[(String, (Constructor[_], Option[AnyRef]))] = {
    // This is the (optional) singleton instance shared by all factory instances
    val singletonInstance: Checked[Option[AnyRef]] = fsConfig.toConfig.getAs[Config]("global")
      .map(c => instantiateSingletonConfig(key, c).toValidated)
      .sequence[ErrorOr, AnyRef]
      .toEither

    (for {
      maybeInstance <- singletonInstance
      parameterTypes = maybeInstance match {
        // If there is a singleton instance, the factory constructor should have a third argument of that type
        case Some(singleton) => List(classOf[Config], classOf[Config], singleton.getClass)
        case None => List(classOf[Config], classOf[Config])
      }
      constructor <- createConstructor(key, fsConfig.toConfig, parameterTypes)
    } yield (key, (constructor, maybeInstance))).toValidated
  }

  // Instantiates the singleton config for a filesystem
  private def instantiateSingletonConfig(filesystem: String, config: Config): Checked[AnyRef] = {
    for {
      constructor <- createConstructor(filesystem, config, List(classOf[Config]))
      instanceConfig = config.getAs[Config]("config").getOrElse(ConfigFactory.empty)
      instance <- Try(constructor.newInstance(instanceConfig)).toChecked
      cast <- instance.cast[AnyRef].toChecked(s"The filesystem global configuration class for $filesystem is not a Java Object")
    } yield cast
  }

  // Create a constructor from a configuration object
  private def createConstructor(key: String, config: Config, parameterTypes: List[Class[_]]): Checked[Constructor[_]] = for {
    clazz <- config.as[Option[String]]("class").toChecked(s"Filesystem configuration $key doesn't have a class field")
    constructor <- createConstructor(key, clazz, parameterTypes)
  } yield constructor

  // Getting a constructor from a class name
  private def createConstructor(filesystem: String, className: String, parameterTypes: List[Class[_]]): Checked[Constructor[_]] = Try (
    Class.forName(className).getConstructor(parameterTypes: _*)
  ).recoverWith({
    case e: ClassNotFoundException => Failure(
      new RuntimeException(s"Class $className for filesystem $filesystem cannot be found in the class path.", e)
    )
    case e: NoSuchMethodException => Failure(
      new RuntimeException(s"Class $className for filesystem $filesystem does not have the required constructor signature: (${parameterTypes.map(_.getCanonicalName).mkString(", ")})", e)
    )
  }).toChecked

  // Instantiate a PathBuilderFactory from its constructor and instance config
  private def instantiate(name: String, constructor: Constructor[_], instanceConfig: Config, global: Option[AnyRef]): Checked[PathBuilderFactory] = {
    for {
      instance <- global match {
        case Some(g) => Try(constructor.newInstance(globalConfig, instanceConfig, g)).toChecked
        case None => Try(constructor.newInstance(globalConfig, instanceConfig)).toChecked
      }
      cast <- instance.cast[PathBuilderFactory].toChecked(s"The filesystem class for $name is not an instance of PathBuilderFactory")
    } yield cast
  }

  // Look for a constructor in the map of known filesystems
  private def getConstructor(fileSystemName: String): Checked[(Constructor[_], Option[AnyRef])] = factoryBuilders
    .get(fileSystemName)
    .toChecked(s"Cannot find a filesystem with name $fileSystemName in the configuration. Available filesystems: ${factoryBuilders.keySet.mkString(", ")}")

  /**
    * Try to find a configured filesystem with the given name and build a PathFactory for it
    * @param name name of the filesystem
    * @param instanceConfig filesystem specific configuration for this instance of the factory to build
    */
  def buildFactory(name: String, instanceConfig: Config): Checked[PathBuilderFactory] = {
    if (DefaultPathBuilderFactory.name.equalsIgnoreCase(name)) DefaultPathBuilderFactory.validNelCheck
    else for {
      constructorAndGlobal <- getConstructor(name)
      factory <- instantiate(name, constructorAndGlobal._1, instanceConfig, constructorAndGlobal._2)
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
      import net.ceedubs.ficus.Ficus._

      def isFilesystemEnabled(configObject: ConfigObject): Boolean = {
        // There may be an "enabled" key on this filesystem. If there is and its value is set to "false" (case insensitive)
        // then ignore this filesystem entry. Otherwise this filesystem is enabled.
        val maybeEnabled = configObject.toConfig.as[Option[String]]("enabled")
        !maybeEnabled.exists(_.equalsIgnoreCase("false"))
      }

      configMap.flatTraverse[ErrorOr, (String, PathBuilderFactory)] {
        // Build the factory for each entry.
        case (key, config: ConfigObject) if isFilesystemEnabled(config) =>
          // This filesystem is enabled, wrap in a List.
          buildFactory(key, config.toConfig).toValidated map { pbf => List(key -> pbf) }
        case (_, _: ConfigObject) =>
          // This filesystem is not enabled, return an empty List.
          List.empty.validNel
        case (key, _) => s"Invalid filesystem backend configuration for $key".invalidNel
      } map { _.toMap } toEither
    } else Map.empty[String, PathBuilderFactory].validNelCheck
  }
}

object CromwellFileSystems {
  val instance = new CromwellFileSystems(ConfigFactory.load())
}
