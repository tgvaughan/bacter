# Dockerfile to build container for unit testing

FROM debian:stable

RUN apt-get update
RUN apt-get install -y openjdk-17-jdk openjfx ant

WORKDIR /root

ADD . ./

ENTRYPOINT JAVA_FX_HOME=/usr/share/java/ ant test
