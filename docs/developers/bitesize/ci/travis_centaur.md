## Travis builds by user access

For infrastructures that require secured credentials, cloud backend tests only run for developers with write access to the broadinstitute/cromwell GitHub. Secure tests are skipped for all other users.

Other backends run tests for any user.

| Backend       | Read-only users | Write/Admin users |
|---------------|:---------------:|:-----------------:|
| AWS           |                 |        ✅         |
| Local         |       ✅        |        ✅         |
| PAPI V2alpha1 |                 |        ✅         |
| PAPI V2beta   |                 |        ✅         |
| SLURM         |       ✅        |        ✅         |
| TES           |       ✅        |        ✅         |

## Upgrade / Horicromtal / etc.

| CI Test Type                  | Cromwell Config                                                  | Centaur Config                                         |
|-------------------------------|------------------------------------------------------------------|--------------------------------------------------------|
| Engine Upgrade                | `(backend)_application.conf`                                     | `centaur_application.conf`*                            |
| Horicromtal                   | `papi_[v2beta or v2alpha1]_horicromtal_application.conf`**       | `centaur_application_`<br>`horicromtal.conf`           |
| Horicromtal<br>Engine Upgrade | `papi_v2beta_application.conf`**                                 | `centaur_application_`<br>`horicromtal_no_assert.conf` |
| Papi Upgrade                  | `papi_v1_v2alpha1_upgrade_application.conf`**                    | `centaur_application.conf`*                            |
| Papi Upgrade<br>New Workflows | `(backend)_application.conf`                                     | `centaur_application.conf`*                            |
| WDL Upgrade                   | `(backend)_application.conf`                                     | `centaur_application.conf`*                            |
| (other)                       | `(backend)_application.conf`                                     | `centaur_application.conf`*                            |

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
\* Centaur Config always uses `centaur_application.conf` except when overridden with `papi_v2alpha1_centaur_application.conf`
or `papi_v2beta_centaur_application.conf`
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
    1. db-mstr: started first
    2. sum-back: runs summarizer
    3. front-back: exposes HTTP
- Horicromtal Engine Upgrade: Combination of Horicromtal and Engine Upgrade
- PAPI Upgrade: Tests run with an older version of Papi and upon restart use a newer version of Papi
- PAPI Upgrade New Workflows: Test definition [does not run any tests](https://travis-ci.org/broadinstitute/cromwell/jobs/475378412)
- WDL Upgrade: Upgrades WDL from draft-2 to 1.0 before testing
- (other): Runs `*.test` files listing the configured backend names

## RDBMS

| Backend | MySQL  | PostgreSQL  | MariaDB  |
|---------|:------:|:-----------:|:--------:|
| AWS     |   ✅   |             |          |
| Local   |   ✅   |      ✅     |          |
| PAPI V2 |   ✅   |             |    ⭕    |
| SLURM   |   ✅   |             |          |
| TES     |   ✅   |             |          |

<small>
⭕ Tests Horicromtal Engine Upgrade versus standard Centaur suite
</small>

All backends run against MySQL. The Local backend also test PostgreSQL, allowing contributors ensure WDLs work with PostgreSQL. MariaDB is tested on a specialized upgrade, where the MySQL connector client is used first, and the MariaDB client is used after restart.
