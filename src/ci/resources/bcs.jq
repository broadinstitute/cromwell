# Convert a bcs date to an jq compatible date.
# The bcs date should already be in UTC.
# https://stedolan.github.io/jq/manual/v1.5/#Dates
def bcsToEpochSeconds(bcs_date):
    bcs_date
    | sub(" "; "T")
    | sub("\\..*$"; "Z")
    | fromdateiso8601;

# Returns true if the jq function `now` is more than `seconds` ahead of `epoch_date_seconds`.
# https://stedolan.github.io/jq/manual/v1.5/#Dates
def isMoreThanSecondsOld(epoch_date_seconds; seconds):
    now - epoch_date_seconds > seconds;

# Filters the bcs date `key` if it is more than `seconds` old.
# https://stedolan.github.io/jq/manual/v1.5/#select(boolean_expression)
def filterMoreThanSecondsOld(key; seconds):
    map(select(isMoreThanSecondsOld(
        bcsToEpochSeconds(key);
        seconds
    )));

# Returns items under `key` that were created more than `seconds` ago.
# For bcs jobs and clusters the key is usually `.`.
# Expects `key` to be able to parse `{ Items: [ .Items[] | {Id, CreationTime} ] }`
def printIdsMoreThanSecondsOld(key; seconds):
    key
    | .Items
    | filterMoreThanSecondsOld(.CreationTime; seconds)
    | .[]
    | .Id;
