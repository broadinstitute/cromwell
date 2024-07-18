For information on Cromwell's Integration Testing Suite, see the [Cromwell documentation on Centaur](https://cromwell.readthedocs.io/en/develop/developers/Centaur/).

### `centaur/src/it`

Classes extending `org.scalatest` that ingest `.test` files and turn them into runnable test suites.

### `centaur/src/main`

#### `/resources`

Collection of `.test` cases. In `test.inc.sh` we map Github Action jobs to case directories with `create_centaur_variables()`. Not all cases are run!

As of July 2024, Centaur searches **recursively** for `.test` files, so they can be placed in subdirectories along with their resources.

#### `/scala`

Functionality to start, stop, and restart the Cromwell server under test. Also contains abstractions for asserting on metadata and workflow outputs.

### `centaur/src/test`

Tests for Centaur itself.
