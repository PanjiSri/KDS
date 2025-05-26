FROM python:3.10.17-slim-bookworm
WORKDIR /usr/src/app
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
COPY ./app ./app
EXPOSE 8050
ENV PYTHONUNBUFFERED=1
CMD ["python", "app/app.py"]