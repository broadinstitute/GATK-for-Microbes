FROM debian:stretch

# apt-get update to fetch available updates
RUN apt-get update -y 

RUN apt-get install default-jdk ant git unzip tabix -y --allow-unauthenticated

RUN git clone https://github.com/RealTimeGenomics/rtg-tools.git && \
    cd rtg-tools && \
    ant runalltests

RUN cd /rtg-tools && \
    ant zip-nojre
RUN unzip /rtg-tools/build/rtg-tools-3.11-39691f9f-base.zip

ENV PATH="$PATH:/rtg-tools-3.11-39691f9f"

