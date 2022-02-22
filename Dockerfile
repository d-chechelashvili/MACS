FROM ubuntu:20.04

RUN apt-get update -y
RUN apt-get install python3.8 -y
RUN apt-get install python3-pip -y
RUN ln -s /usr/bin/python3 /usr/bin/python

WORKDIR /sandbox

CMD pip3 install -r requirements.txt & /bin/bash