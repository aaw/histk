FROM alpine

RUN apk add --no-cache bash build-base git linux-headers ruby ruby-bundler ruby-io-console

RUN git clone https://github.com/antirez/redis.git /opt/redis && \
    cd /opt/redis && \
    git checkout unstable && \
    make && \
    cd -

RUN gem install --no-ri --no-rdoc test-unit redis table_print

COPY src test data /opt/histk/

WORKDIR /opt/histk

RUN make clean && make

ENTRYPOINT ["/bin/bash", "/opt/histk/entrypoint.sh"]

CMD "/bin/bash"