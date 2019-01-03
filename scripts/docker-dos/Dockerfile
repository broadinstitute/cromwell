FROM openjdk:8

RUN mkdir ~/.ammonite && curl -Ls -o ~/.ammonite/predef.scala https://git.io/vro0a
RUN curl -Ls -o /bin/amm https://git.io/vdNv2 && chmod +x /bin/amm
RUN mkdir /scripts

ADD dosUrlLocalizer.sc /scripts
