
# Use an official Python runtime as a parent image
#FROM python:3.9

# Set the working directory in the container
#WORKDIR /app

# Copy the current directory contents into the container at /app
#COPY . /app

#RUN pip3 install scikit-image
#RUN pip3 install readlif
#RUN pip3 install matplotlib

# Command to run the Python script
#CMD ["python", "script_docker.py", "first_frame", "/app/movies.lif"]

FROM ubuntu:20.04

RUN apt-get update
RUN apt-get install -y gcc make libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev 
RUN apt-get install -y bzip2

#python3 installation
RUN apt-get install build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev -y
COPY ./Python-3.11.2.tar.xz /tmp/Python-3.11.2.tar.xz
RUN tar -xf /tmp/Python-3.11.2.tar.xz

WORKDIR /Python-3.11.2
RUN ./configure --enable-optimizations
RUN make install

RUN mkdir /app; chmod 777 /app


# Copy the current directory contents into the container at /app

COPY . /app

RUN pip3 install scikit-image
RUN pip3 install readlif
RUN pip3 install matplotlib


WORKDIR /app

CMD ["bash"]