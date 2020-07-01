FROM centos:7

RUN yum update -y && yum install -y cmake wget gcc gcc-c++ bzip2 autoconf zlib-devel make zip unzip curl-devel

ADD entry.sh /entry.sh
RUN chmod a+x /entry.sh

ENTRYPOINT ["/entry.sh"]
