FROM python:3.10.17-slim-bookworm
WORKDIR /usr/src/app
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
COPY ./app ./app
ENV PYTHONPATH=/usr/src/app:$PYTHONPATH
EXPOSE 8050
ENV PYTHONUNBUFFERED=1
ENV PYTHONWARNINGS="ignore::UserWarning"
CMD ["python", "-m", "app.app"]