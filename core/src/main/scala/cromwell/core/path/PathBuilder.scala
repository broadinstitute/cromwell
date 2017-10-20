package cromwell.core.path

import scala.util.Try

trait PathBuilder {
  def name: String

  def build(pathAsString: String): Try[Path]
}

/**
  * A path that was built by a PathBuilder.
  *
  * Encapsulates an underlying java nio path adding utility methods.
  *
  * @see [[cromwell.core.path.NioPathMethods]]
  * @see [[cromwell.core.path.BetterFileMethods]]
  * @see [[cromwell.core.path.EvenBetterPathMethods]]
  */
trait Path extends PathObjectMethods with NioPathMethods with BetterFileMethods with EvenBetterPathMethods {
  /**
    * A reference to the underlying nioPath, used to create new java.nio.Path's that will then be sent to newPath
    * for wrapping.
    *
    * @return The internal nioPath.
    */
  protected def nioPath: NioPath

  /**
    * Wraps the nioPath in our container type.
    *
    * The nioPath should be from the same filesystem.
    *
    * @param nioPath The nioPath to be wrapped.
    * @return A new Path.
    */
  protected def newPath(nioPath: NioPath): Path

  /**
    * Returns the path as a string. This path must be usable by this path builder or another in the list to build paths,
    * and when isAbsolute returns true, should be usable by external users for accessing their files.
    *
    * Examples:
    * {{{
    * path                          path.pathAsString
    * ---------------------------   ---------------------------
    * gs://bucket/path/to/file      gs://bucket/path/to/file
    * gs://bucket/path/to/my file   gs://bucket/path/to/my file
    * /mount/path/to/file           /mount/path/to/file
    * /mount/path/to/my file        /mount/path/to/my file
    * }}}
    *
    * If isAbsolute returns false, the behavior may be different depending on the Path/PathBuilder:
    * {{{
    * path                      path.isAbsolute  path.getName(2)      p.gN(2).isAbsolute     p.gN(2).pathAsString
    * ------------------------  ---------------  ------------------   ---------------------  -------------------------
    * gs://bucket/path/to/file  true             gs://bucket/file     false                  gs://bucket/file
    * /mount/path/to/file       true             /mount/path/to/file  false                  to
    * mount/path/to/file        false            mount/path/to/file   false                  to
    * }}}
    *
    * @return Path as a string
    */
  def pathAsString: String

  /**
    * Returns the unencoded uri of the string of the string stripping off the scheme.
    *
    * This represents the old model of {{{path.toUri.getHost + path.toUri.getPath}}}, but doesn't encode the string.
    *
    * Examples:
    * {{{
    * path                          path.pathWithoutScheme
    * ---------------------------   -----------------------
    * gs://bucket/path/to/file      bucket/path/to/file
    * gs://bucket/path/to/my file   bucket/path/to/my file
    * /mount/path/to/file           /mount/path/to/file
    * /mount/path/to/my file        /mount/path/to/my file
    * }}}
    *
    * This path is not relative to the root of the underlying path:
    * {{{
    * path                       path.root      path.pathWithoutScheme   path.root.resolve(path.pathWithoutScheme)
    * ------------------------   ------------   ----------------------   -------------------------------------------
    * gs://bucket/path/to/file   gs://bucket/   bucket/path/to/file      gs://bucket/bucket/path/to/file  (BAD!!)
    * /mount/path/to/file        /              /mount/path/to/file      /mount/path/to/file  (ok, but still don't.)
    * }}}
    *
    * The returned value is also different from a [[java.net.URI]], in that the value should not be encoded:
    * {{{
    * path                          path.pathWithoutScheme    path.toUri
    * ---------------------------   -----------------------   ------------------------------------------
    * gs://bucket/path/to/file      bucket/path/to/file       gs://bucket/path/to/file
    * gs://bucket/path/to/my file   bucket/path/to/my file    gs://bucket/path/to/my%20file  (BAD!!)
    * /mount/path/to/file           /mount/path/to/file       file:///mount/path/to/file
    * /mount/path/to/my file        /mount/path/to/my file    file:///mount/path/to/my%20file  (BAD!!)
    * }}}
    *
    * If isAbsolute returns false, the behavior may be different depending on the Path/PathBuilder:
    * {{{
    * path                      path.isAbsolute  path.getName(2)   p.gN(2).isAbsolute     p.gN(2).pathWithoutScheme
    * ------------------------  ---------------  ----------------  ---------------------  -------------------------
    * gs://bucket/path/to/file  true             gs://bucket/file  false                  bucket/file
    * /mount/path/to/file       true             to                false                  /<pwd>/to
    * }}}

    *
    * @return A "relative" path
    */
  def pathWithoutScheme: String

  // Used by various extension traits within this scala package
  private[path] final def nioPathPrivate: NioPath = nioPath

  // Used within BetterFileMethods
  private[path] final def betterFile: better.files.File = nioPathPrivate

  // Some Path methods return null.
  private[path] final def newPathOrNull(nioPath: NioPath) = Option(nioPath).map(newPath).orNull
}
