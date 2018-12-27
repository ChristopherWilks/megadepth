FROM gcc:5.5

RUN apt-get update -y && apt-get install -y cmake zip unzip

RUN wget http://packages.seqan.de/seqan-library/seqan-library-2.4.0.zip -O /tmp/seqan-library-2.4.0.zip
RUN cd /tmp && unzip seqan-library-2.4.0.zip && \
    cp -r seqan-library-2.4.0/include/* /usr/local/include/ && \
    cp -r seqan-library-2.4.0/share/* /usr/local/share/ && \
    rm -rf seqan-library-2.4.0

ADD entry.sh /entry.sh
RUN chmod a+x /entry.sh

ENTRYPOINT ["/entry.sh"]
