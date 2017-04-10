#!/bin/sh

IMAGE_NAME=shu65/astralroc
IMAGE=${IMAGE_NAME}:latest

echo "image = ${IMAGE}"
cd ../
docker build -t ${IMAGE} -f ./docker/Dockerfile .
