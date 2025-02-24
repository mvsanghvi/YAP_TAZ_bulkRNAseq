FROM python:3

# Get latest STAR source from releases
RUN git clone https://github.com/alexdobin/STAR.git

RUN apt-get update && \
    apt-get install g++ && \
    apt-get install make && \
    apt-get install xxd

# Compile STAR
RUN cd STAR/source && \
    make clean && \
    make STAR

WORKDIR /STAR/source

CMD ["/bin/bash"]
