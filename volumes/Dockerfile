FROM python:3.9.5-slim
RUN mkdir /app
WORKDIR /app
ADD requirements.txt ./requirements.txt 

RUN apt-get clean && apt-get update \
   && apt-get install -y  --no-install-recommends \
   g++ \
   gcc \
   gfortran \
   libgfortran5 \
   rsync \
   && apt-get install -y make build-essential

RUN pip install --default-timeout=150 -r requirements.txt  


RUN  ln -snf /usr/share/zoneinfo/America/Lima /etc/localtime \
   && echo "America/Lima" > /etc/timezone \
   && dpkg-reconfigure -f noninteractive tzdata \ 
   && rm -rf /var/lib/apt/lists/*