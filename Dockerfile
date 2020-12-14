FROM python:latest

# apt-get update to fetch available updates
RUN pip3 install numpy
RUN pip3 install CrossMap && pip3 install CrossMap --upgrade

