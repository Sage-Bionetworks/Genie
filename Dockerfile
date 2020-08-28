FROM python:3.7

COPY . .
RUN python setup.py install
