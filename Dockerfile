# This Dockerfile builds a Python dev environment
FROM ubuntu

RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-venv \
    build-essential \
    && pip install twine build

COPY . /code
WORKDIR /code
RUN python3 -m build

#CMD ["python3", "-m", "twine", "upload dist/*"]
