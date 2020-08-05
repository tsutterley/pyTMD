FROM python:3.7-alpine

LABEL "com.github.actions.name"="pyTMD"
LABEL "com.github.actions.description"="Run pytest commands"
LABEL "com.github.actions.icon"="upload-cloud"
LABEL "com.github.actions.color"="yellow"

RUN apk add --no-cache bash
RUN pip install --upgrade pip
RUN pip install pytest pytest-cov
RUN pip install -r requirements.txt
RUN python --version ; pip --version ; pytest --version
