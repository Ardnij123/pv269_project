FROM python:3.13.3-slim

RUN apt-get update
RUN apt-get install -y g++ libz-dev make libcairo2

RUN python -m pip install --upgrade pip

ADD https://github.com/marbl/ModDotPlot.git /usr/src/moddotplot
RUN cd /usr/src/moddotplot && python -m pip install . && rm -rf /root/.cache

RUN apt-get autoremove

CMD ["moddotplot"]
