workflow Hello {

    call HelloWorld

    output {
        String message=HelloWorld.out
    }
}

task HelloWorld{
    command {
        echo "hello"
    }
    output {
        String out = stdout()
    }
}

task HelloWorld{
    command {
        echo "goodbye"
    }
    output {
        String out = stdout()
    }
}
