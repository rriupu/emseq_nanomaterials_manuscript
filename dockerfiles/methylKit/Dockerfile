FROM r-base:4.5.1

COPY requirements-deb.txt /opt/
COPY requirements.R /opt/
WORKDIR /opt/

RUN apt-get update -y
RUN xargs apt-get install -y < requirements-deb.txt \
    && Rscript requirements.R