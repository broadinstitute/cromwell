##
## A set of expressions which would cause errors, except luckily they're short-circuited
##
workflow short_circuit {

  if (false && 1/0 == 0) {
    Int a = 1
  }

  if (true || 1/0 == 0) {
    Int b = 1
  }

  Int c = if true then 1 else 1/0
  Int d = if false then 1/0 else 1

}
