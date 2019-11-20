## Travis builds by user access

For infrastructures that require secured credentials, cloud backend tests only run for developers with write access to the broadinstitute/cromwell GitHub. Secure tests are skipped for all other users.

Other backends run tests for any user.

| Backend | Read-only users | Write/Admin users |
|---------|:---------------:|:-----------------:|
| AWS     |                 |        ✅         |
| BCS     |                 |        ✅         |
| Local   |       ✅        |        ✅         |
| PAPI V1 |                 |        ✅         |
| PAPI V2 |                 |        ✅         |
| SLURM   |       ✅        |        ✅         |
| TES     |       ✅        |        ✅         |

## Upgrade / Horicromtal / etc.

| CI Test Type                  | Cromwell Config                          | Centaur Config                                         |
|-------------------------------|------------------------------------------|--------------------------------------------------------|
| Engine Upgrade                | `(backend)_application.conf`             | `centaur_application.conf`*                            |
| Horicromtal                   | `papi_v2_horicromtal_application.conf`** | `centaur_application_`<br>`horicromtal.conf`           |
| Horicromtal<br>Engine Upgrade | `papi_v2_application.conf`**             | `centaur_application_`<br>`horicromtal_no_assert.conf` |
| Papi Upgrade                  | `papi_v1_v2_upgrade_application.conf`**  | `centaur_application.conf`*                            |
| Papi Upgrade<br>New Workflows | `(backend)_application.conf`             | `centaur_application.conf`*                            |
| WDL Upgrade                   | `(backend)_application.conf`             | `centaur_application.conf`*                            |
| (other)                       | `(backend)_application.conf`             | `centaur_application.conf`*                            |

| CI Test Type                  | ScalaTest Spec              | Test Directory                      |
|-------------------------------|-----------------------------|-------------------------------------|
| Engine Upgrade                | `EngineUpgradeTestCaseSpec` | `engineUpgradeTestCases`            |
| Horicromtal                   | `CentaurTestSuite`          | `standardTestCases`***              |
| Horicromtal<br>Engine Upgrade | `EngineUpgradeTestCaseSpec` | `engineUpgradeTestCases`***         |
| PAPI Upgrade                  | `PapiUpgradeTestCaseSpec`   | `papiUpgradeTestCases`              |
| PAPI Upgrade<br>New Workflows | `CentaurTestSuite`          | `papiUpgradeNewWorkflowsTestCases`  |
| WDL Upgrade                   | `WdlUpgradeTestCaseSpec`    | `standardTestCases`****             |
| (other)                       | `CentaurTestSuite`          | `standardTestCases`                 |

<small>
\* Centaur Config always uses `centaur_application.conf` except when overridden with `papi_v2_centaur_application.conf`
  ([48 preview link](https://github.com/broadinstitute/cromwell/blob/a7d0601/src/ci/bin/test.inc.sh#L455-L457))  
\*\* Cromwell Config overrides
  ([47 link](https://github.com/broadinstitute/cromwell/blob/47/src/ci/bin/test.inc.sh#L213-L221))  
\*\*\* Test Directory overrides
  ([47 link](https://github.com/broadinstitute/cromwell/blob/47/src/ci/bin/test.inc.sh#L440-L449))  
\*\*\*\* Test Directory only tests tagged with `wdl_upgrade`
  ([47 link](https://github.com/broadinstitute/cromwell/blob/47/centaur/src/main/resources/standardTestCases/write_lines.test#L3))  
</small>

- Engine Upgrade: Retrieves the [Cromwell Version](https://github.com/broadinstitute/cromwell/blob/47/project/Version.scala#L8) then retrieves the previous jar/docker-image from DockerHub. Centaur starts with the prior version, then restarts with the compiled source code.
- Horicromtal: Runs a [docker-compose](https://github.com/broadinstitute/cromwell/blob/47/src/ci/docker-compose/docker-compose-horicromtal.yml) with:
    1. cromwell-database-master: started first
    2. cromwell-summarizer-plus-backend: runs summarizer
    3. cromwell-frontend-plus-backend: exposes HTTP
- Horicromtal Engine Upgrade: Combination of Horicromtal and Engine Upgrade
- PAPI Upgrade: Tests run with Papi V1 and upon restart use Papi V2
- PAPI Upgrade New Workflows: Test definition [does not run any tests](https://travis-ci.org/broadinstitute/cromwell/jobs/475378412)
- WDL Upgrade: Upgrades WDL from draft-2 to 1.0 before testing
- (other): Runs `*.test` files listing the configured backend names

## RDBMS

| Backend | MySQL  | PostgreSQL  | MariaDB  |
|---------|:------:|:-----------:|:--------:|
| AWS     |   ✅   |             |          |
| BCS     |   ✅   |             |          |
| Local   |   ✅   |      ✅     |          |
| PAPI V1 |   ✅   |             |          |
| PAPI V2 |   ✅   |             |    ⭕    |
| SLURM   |   ✅   |             |          |
| TES     |   ✅   |             |          |

<small>
⭕ Tests Horicromtal Engine Upgrade versus standard Centaur suite
</small>

All backends run against MySQL. The Local backend also test PostgreSQL, allowing contributors ensure WDLs work with PostgreSQL. MariaDB is tested on a specialized upgrade, where the MySQL connector client is used first, and the MariaDB client is used after restart.
