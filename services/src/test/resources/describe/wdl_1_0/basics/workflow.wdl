version 1.0

workflow basics {

  meta {
    email: "skroob@spaceballs.gov"
    author: "President Skroob"
    description: "Spaceballs: The Unit Test"
  }

  input {
    Int i = 2 + 2
    Int j = 4
    File f
  }
}
