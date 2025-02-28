version 1.0

task makeAnEmoji {
    command {
        AN_EMOJI='\360\237\246\204\n'
        echo "This is my stdout file"
        echo "This is my stderr file" >&2
        printf $AN_EMOJI >&2
    }
    runtime {
        docker: 'ubuntu:latest'
    }
}

workflow emoji_in_stderr {
    call makeAnEmoji {}
}