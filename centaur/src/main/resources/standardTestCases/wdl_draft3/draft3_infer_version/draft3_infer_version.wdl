
# Some heading commentary...

# FIXME/TODO:
#
# WDL version metadata is casted to String somewhere between WDL parsing and metadata checking in Centaur.
# As a temporary measure to make this test pass, this test is expecting the WDL version to be a String.
# This is an issue in Centaur and does not affect production.
version 1.0

workflow draft3_infer_version {
  input {
    Int i = 12
  }

  output {
    Int j = i * 55
  }
}
