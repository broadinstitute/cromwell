[![Build Status](https://travis-ci.org/broadinstitute/lenthall.svg?branch=develop)](https://travis-ci.org/broadinstitute/lenthall?branch=develop)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/lenthall/badge.svg?branch=develop)](https://coveralls.io/r/broadinstitute/lenthall?branch=develop)

Lenthall
========

[Cromwell](https://github.com/broadinstitute/cromwell) Common Code

# Mailing List

The [Cromwell Mailing List](https://groups.google.com/a/broadinstitute.org/forum/?hl=en#!forum/cromwell) is cromwell@broadinstitute.org.

If you have any questions, suggestions or support issues please send them to this list. To subscribe you can either join via the link above or send an email to cromwell+subscribe@broadinstitute.org.

# Requirements

The following is the toolchain used for development of Lenthall.  Other versions may work, but these are recommended.

* [Scala 2.11.7](http://www.scala-lang.org/news/2.11.7)
* [SBT 0.13.8](https://github.com/sbt/sbt/releases/tag/v0.13.8)
* [Java 8](http://www.oracle.com/technetwork/java/javase/overview/java8-2100321.html)

# Building

`sbt compile` will build a library JAR in `target/scala-2.11/`

Tests are run via `sbt test`.
