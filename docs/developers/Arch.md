# Database Reference


## All tables

|Table|Description|
|--|--|
|Call Caching Aggregation| |
|Call Caching Detritus| |
|Call Caching| |
|Call Caching Hash| |
|Call Caching Simpleton| |
|Custom Label| |
|Docker Hash Store| |
|Job Key Value| |
|Job Store| |
|Job Store Simpleton| |
|Metadata Entry| |
|Sub Workflow Store| |
|Summary Status| |
|Workflow Metadata Summary| |
|Workflow Store| |

## Call Caching Aggregation

| Field                             | Type         | Null | Key | Default | Extra          |
|-----------------------------------|--------------|------|-----|---------|----------------|
| CALL_CACHING_AGGREGATION_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| CALL_CACHING_ENTRY_ID             | int(11)      | NO   | MUL | NULL    |                |
| BASE_AGGREGATION                  | varchar(255) | NO   | MUL | NULL    |                |
| INPUT_FILES_AGGREGATION           | varchar(255) | YES  |     | NULL    |                |

## Call Caching Detritus

| Field                          | Type         | Null | Key | Default | Extra          |
|--------------------------------|--------------|------|-----|---------|----------------|
| CALL_CACHING_DETRITUS_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| DETRITUS_KEY                   | varchar(255) | YES  |     | NULL    |                |
| DETRITUS_VALUE                 | longtext     | YES  |     | NULL    |                |
| CALL_CACHING_ENTRY_ID          | int(11)      | YES  | MUL | NULL    |                |

## Call Caching Entry

| Field                     | Type         | Null | Key | Default | Extra          |
|---------------------------|--------------|------|-----|---------|----------------|
| CALL_CACHING_ENTRY_ID     | int(11)      | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID   | varchar(255) | YES  | MUL | NULL    |                |
| CALL_FULLY_QUALIFIED_NAME | varchar(255) | YES  |     | NULL    |                |
| JOB_INDEX                 | int(11)      | YES  |     | NULL    |                |
| RETURN_CODE               | int(11)      | YES  |     | NULL    |                |
| ALLOW_RESULT_REUSE        | tinyint(4)   | YES  |     | 1       |                |
| JOB_ATTEMPT               | int(11)      | YES  |     | NULL    |                |

## Call Caching Hash Entry

| Field                      | Type         | Null | Key | Default | Extra          |
|----------------------------|--------------|------|-----|---------|----------------|
| CALL_CACHING_HASH_ENTRY_ID | bigint(20)   | NO   | PRI | NULL    | auto_increment |
| HASH_KEY                   | varchar(255) | NO   |     | NULL    |                |
| HASH_VALUE                 | varchar(255) | NO   |     | NULL    |                |
| CALL_CACHING_ENTRY_ID      | int(11)      | YES  | MUL | NULL    |                |

## Call Caching Simpleton Entry

| Field                           | Type         | Null | Key | Default | Extra          |
|---------------------------------|--------------|------|-----|---------|----------------|
| CALL_CACHING_SIMPLETON_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| HASH_KEY                        | varchar(255) | NO   |     | NULL    |                |
| HASH_VALUE                      | varchar(255) | NO   |     | NULL    |                |
| CALL_CACHING_ENTRY_ID           | int(11)      | YES  | MUL | NULL    |                |

## Custom Labels

| Field                   | Type         | Null | Key | Default | Extra          |
|-------------------------|--------------|------|-----|---------|----------------|
| CUSTOM_LABEL_ENTRY_ID   | bigint(20)   | NO   | PRI | NULL    | auto_increment |
| CUSTOM_LABEL_KEY        | varchar(255) | YES  | MUL | NULL    |                |
| CUSTOM_LABEL_VALUE      | varchar(255) | YES  |     | NULL    |                |
| WORKFLOW_EXECUTION_UUID | varchar(100) | NO   | MUL | NULL    |                |

## Docker hash store

| Field                      | Type         | Null | Key | Default | Extra          |
|----------------------------|--------------|------|-----|---------|----------------|
| DOCKER_HASH_STORE_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID    | varchar(255) | NO   | MUL | NULL    |                |
| DOCKER_TAG                 | varchar(255) | NO   |     | NULL    |                |
| DOCKER_HASH                | varchar(255) | NO   |     | NULL    |                |


## Job Key Value

| Field                     | Type         | Null | Key | Default | Extra          |
|---------------------------|--------------|------|-----|---------|----------------|
| JOB_KEY_VALUE_ENTRY_ID    | int(11)      | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID   | varchar(255) | NO   | MUL | NULL    |                |
| CALL_FULLY_QUALIFIED_NAME | varchar(255) | YES  |     | NULL    |                |
| JOB_INDEX                 | int(11)      | YES  |     | NULL    |                |
| JOB_ATTEMPT               | int(11)      | YES  |     | NULL    |                |
| STORE_KEY                 | varchar(255) | NO   |     | NULL    |                |
| STORE_VALUE               | varchar(255) | NO   |     | NULL    |                |

## Sub Workflow Store

| Field                          | Type         | Null | Key | Default | Extra          |
|--------------------------------|--------------|------|-----|---------|----------------|
| SUB_WORKFLOW_STORE_ENTRY_ID    | int(11)      | NO   | PRI | NULL    | auto_increment |
| ROOT_WORKFLOW_ID               | int(11)      | NO   | MUL | NULL    |                |
| PARENT_WORKFLOW_EXECUTION_UUID | varchar(255) | NO   | MUL | NULL    |                |
| CALL_FULLY_QUALIFIED_NAME      | varchar(255) | NO   |     | NULL    |                |
| CALL_INDEX                     | int(11)      | NO   |     | NULL    |                |
| CALL_ATTEMPT                   | int(11)      | NO   |     | NULL    |                |
| SUB_WORKFLOW_EXECUTION_UUID    | varchar(255) | NO   |     | NULL    |                |

## Summary Status

| Field                   | Type         | Null | Key | Default | Extra          |
|-------------------------|--------------|------|-----|---------|----------------|
| SUMMARY_STATUS_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| SUMMARY_TABLE_NAME      | varchar(255) | NO   | MUL | NULL    |                |
| SUMMARIZED_TABLE_NAME   | varchar(255) | NO   |     | NULL    |                |
| MAXIMUM_ID              | bigint(20)   | NO   |     | NULL    |                |

## Job Store

| Field                     | Type         | Null | Key | Default | Extra          |
|---------------------------|--------------|------|-----|---------|----------------|
| JOB_STORE_ENTRY_ID        | int(11)      | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID   | varchar(255) | YES  | MUL | NULL    |                |
| CALL_FULLY_QUALIFIED_NAME | varchar(255) | YES  |     | NULL    |                |
| JOB_INDEX                 | int(11)      | YES  |     | NULL    |                |
| JOB_ATTEMPT               | int(11)      | YES  |     | NULL    |                |
| JOB_SUCCESSFUL            | tinyint(4)   | NO   |     | NULL    |                |
| RETURN_CODE               | int(11)      | YES  |     | NULL    |                |
| EXCEPTION_MESSAGE         | longtext     | YES  |     | NULL    |                |
| RETRYABLE_FAILURE         | tinyint(4)   | YES  |     | NULL    |                |


## Job Store Simpleton

| Field                        | Type         | Null | Key | Default | Extra          |
|------------------------------|--------------|------|-----|---------|----------------|
| JOB_STORE_SIMPLETON_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| SIMPLETON_KEY                | varchar(255) | NO   |     | NULL    |                |
| SIMPLETON_VALUE              | longtext     | YES  |     | NULL    |                |
| WDL_TYPE                     | varchar(255) | NO   |     | NULL    |                |
| JOB_STORE_ENTRY_ID           | int(11)      | YES  | MUL | NULL    |                |


## Metadata

| Field                   | Type         | Null | Key | Default | Extra          |
|-------------------------|--------------|------|-----|---------|----------------|
| METADATA_JOURNAL_ID     | bigint(20)   | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID | varchar(255) | NO   | MUL | NULL    |                |
| METADATA_KEY            | varchar(255) | NO   |     | NULL    |                |
| CALL_FQN                | varchar(255) | YES  |     | NULL    |                |
| JOB_SCATTER_INDEX       | int(11)      | YES  |     | NULL    |                |
| JOB_RETRY_ATTEMPT       | int(11)      | YES  |     | NULL    |                |
| METADATA_VALUE          | longtext     | YES  |     | NULL    |                |
| METADATA_TIMESTAMP      | datetime     | NO   |     | NULL    |                |
| METADATA_VALUE_TYPE     | varchar(10)  | YES  |     | NULL    |                |


## Workflow Metadata Summary

| Field                              | Type         | Null | Key | Default | Extra          |
|------------------------------------|--------------|------|-----|---------|----------------|
| WORKFLOW_METADATA_SUMMARY_ENTRY_ID | bigint(20)   | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID            | varchar(100) | NO   | UNI | NULL    |                |
| WORKFLOW_NAME                      | varchar(100) | YES  | MUL | NULL    |                |
| WORKFLOW_STATUS                    | varchar(50)  | YES  | MUL | NULL    |                |
| START_TIMESTAMP                    | datetime     | YES  |     | NULL    |                |
| END_TIMESTAMP                      | datetime     | YES  |     | NULL    |                |
| SUBMISSION_TIMESTAMP               | datetime     | YES  |     | NULL    |                |

## Workflow Store

| Field                   | Type         | Null | Key | Default | Extra          |
|-------------------------|--------------|------|-----|---------|----------------|
| WORKFLOW_STORE_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID | varchar(255) | YES  | UNI | NULL    |                |
| WORKFLOW_DEFINITION     | longtext     | YES  |     | NULL    |                |
| WORKFLOW_INPUTS         | longtext     | YES  |     | NULL    |                |
| WORKFLOW_OPTIONS        | longtext     | YES  |     | NULL    |                |
| WORKFLOW_STATE          | varchar(20)  | YES  | MUL | NULL    |                |
| SUBMISSION_TIME         | datetime     | NO   |     | NULL    |                |
| IMPORTS_ZIP             | longblob     | YES  |     | NULL    |                |
| CUSTOM_LABELS           | longtext     | NO   |     | NULL    |                |
| WORKFLOW_TYPE           | varchar(30)  | YES  |     | NULL    |                |
| WORKFLOW_TYPE_VERSION   | varchar(255) | YES  |     | NULL    |                |
| WORKFLOW_ROOT           | varchar(100) | YES  |     | NULL    |                |
| CROMWELL_ID             | varchar(100) | YES  |     | NULL    |                |
| HEARTBEAT_TIMESTAMP     | timestamp    | YES  |     | NULL    |                |

