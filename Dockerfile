FROM ubuntu:20.04

RUN apt update
RUN apt upgrade

# Install python 2.7
RUN apt install curl -y
RUN apt install python2.7 -y

# Install python2 pip
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
RUN python2.7 get-pip.py

# Alias python to python2.7
RUN ln -s /usr/bin/python2.7 /usr/bin/python

WORKDIR /sandbox

CMD /bin/bash
