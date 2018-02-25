task put_that_coffee_down {
  command {
    echo "coffee's for valid options json only"
  }
}

workflow invalid_options_json {
  call put_that_coffee_down
}
